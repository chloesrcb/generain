#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.spatial.distance import pdist, squareform
from scipy.stats import rankdata, norm
from numba import njit


#%%
###############################################################################
# 1. DISTANCE MATRIX
###############################################################################

def get_dist_mat_from_lonlat(df_loc, lon_col="Longitude", lat_col="Latitude"):
    """
    Approximation locale en mètres à partir des lon/lat.
    Si tu veux être parfaitement cohérente avec Lambert-93, on pourra remplacer.
    """
    lon = np.asarray(df_loc[lon_col], dtype=float)
    lat = np.asarray(df_loc[lat_col], dtype=float)

    lat0 = np.mean(lat)
    x = lon * 111320.0 * np.cos(np.deg2rad(lat0))
    y = lat * 110540.0

    coords = np.column_stack([x, y])
    dist_mat = squareform(pdist(coords))
    return dist_mat


#%%
###############################################################################
# 2. HELPERS
###############################################################################

def make_radius_quantile(dist_mat, n_bins=10):
    d = dist_mat[np.triu_indices_from(dist_mat, k=1)]
    breaks = np.quantile(d, np.linspace(0, 1, n_bins + 1))
    breaks = np.unique(breaks)
    return breaks


def get_midpoints(breaks):
    return 0.5 * (breaks[:-1] + breaks[1:])


def vario_spatemp(chi):
    chi = np.clip(chi, 1e-6, 1 - 1e-6)
    return 2.0 * norm.ppf(1 - 0.5 * chi) ** 2


def check_station_alignment(location_gauges, rain_df, station_col="codestation"):
    loc_names = list(location_gauges[station_col].astype(str))
    rain_names = list(rain_df.columns.astype(str))

    print("Nb stations location_gauges:", len(loc_names))
    print("Nb colonnes pluie:", len(rain_names))

    only_loc = [x for x in loc_names if x not in rain_names]
    only_rain = [x for x in rain_names if x not in loc_names]

    if only_loc:
        print("Stations dans location_gauges mais pas dans rain:", only_loc)
    if only_rain:
        print("Stations dans rain mais pas dans location_gauges:", only_rain)

    print("Ordre identique :", loc_names == rain_names)


#%%
###############################################################################
# 3. GLOBAL UNIFORM MARGINS
###############################################################################

def make_global_uniform_margins(data_rain):
    """
    Transforme chaque station séparément en pseudo-uniforme global :
        U_{s,t} = rank(X_{s,t}) / (n_s + 1)
    où n_s = nb d'observations non manquantes de la station s.

    Important :
    - le ranking est fait une seule fois par station
    - pas localement par paire/lag
    """
    if isinstance(data_rain, pd.DataFrame):
        data_raw = data_rain.to_numpy(dtype=np.float64)
        colnames = list(data_rain.columns)
    else:
        data_raw = np.asarray(data_rain, dtype=np.float64)
        colnames = [f"s{i}" for i in range(data_raw.shape[1])]

    n_time, n_sites = data_raw.shape
    data_unif = np.full((n_time, n_sites), np.nan, dtype=np.float64)

    for j in range(n_sites):
        x = data_raw[:, j]
        ok = ~np.isnan(x)
        x_ok = x[ok]

        if x_ok.size == 0:
            continue

        ranks = rankdata(x_ok, method="average")
        data_unif[ok, j] = ranks / (x_ok.size + 1.0)

    return data_raw, data_unif, colnames


#%%
###############################################################################
# 4. FAST CHI COMPUTATION WITH GLOBAL MARGINS
###############################################################################

@njit
def _compute_chi_for_pair_tau_global_u(x1_raw, x2_raw, u1_full, u2_full, tau, q, remove_zeros):
    """
    Calcule chi pour une paire (s1, s2) et un lag tau :
      chi = P(U2 > q | U1 > q)
    avec :
      - U1, U2 = marges uniformes globales
      - filtrage après alignement temporel
      - retrait optionnel des couples (0,0) sur données brutes
    """
    n_time = len(x1_raw)
    m = n_time - tau

    if m < 2:
        return np.nan, 0

    count = 0
    for t in range(m):
        a_raw = x1_raw[t]
        b_raw = x2_raw[t + tau]
        a_u = u1_full[t]
        b_u = u2_full[t + tau]

        if np.isnan(a_raw) or np.isnan(b_raw) or np.isnan(a_u) or np.isnan(b_u):
            continue

        if remove_zeros and (a_raw == 0.0 and b_raw == 0.0):
            continue

        count += 1

    if count < 2:
        return np.nan, count

    n_exc = 0
    n_joint = 0

    for t in range(m):
        a_raw = x1_raw[t]
        b_raw = x2_raw[t + tau]
        a_u = u1_full[t]
        b_u = u2_full[t + tau]

        if np.isnan(a_raw) or np.isnan(b_raw) or np.isnan(a_u) or np.isnan(b_u):
            continue

        if remove_zeros and (a_raw == 0.0 and b_raw == 0.0):
            continue

        if a_u > q:
            n_exc += 1
            if b_u > q:
                n_joint += 1

    if n_exc == 0:
        return np.nan, count

    return n_joint / n_exc, count


@njit
def _calc_chi_spatemp_global_u_numba(data_raw, data_unif, pair_i, pair_j, pair_dist, tau_max, q, remove_zeros):
    n_pairs = len(pair_i)
    n_rows = n_pairs * (tau_max + 1)

    s1_out = np.empty(n_rows, dtype=np.int64)
    s2_out = np.empty(n_rows, dtype=np.int64)
    tau_out = np.empty(n_rows, dtype=np.int64)
    dist_out = np.empty(n_rows, dtype=np.float64)
    chi_out = np.empty(n_rows, dtype=np.float64)
    n_used_out = np.empty(n_rows, dtype=np.int64)

    row = 0
    for p in range(n_pairs):
        i = pair_i[p]
        j = pair_j[p]
        d = pair_dist[p]

        x1_raw = data_raw[:, i]
        x2_raw = data_raw[:, j]
        u1_full = data_unif[:, i]
        u2_full = data_unif[:, j]

        for tau in range(tau_max + 1):
            chi_val, n_used = _compute_chi_for_pair_tau_global_u(
                x1_raw, x2_raw, u1_full, u2_full, tau, q, remove_zeros
            )

            s1_out[row] = i
            s2_out[row] = j
            tau_out[row] = tau
            dist_out[row] = d
            chi_out[row] = chi_val
            n_used_out[row] = n_used
            row += 1

    return s1_out, s2_out, tau_out, dist_out, chi_out, n_used_out


def calc_chi_spatemp_fast_global_u(data_rain, dist_mat, tau_max, q, remove_zeros=True):
    data_raw, data_unif, colnames = make_global_uniform_margins(data_rain)

    iu = np.triu_indices(dist_mat.shape[0], k=1)
    pair_i = iu[0].astype(np.int64)
    pair_j = iu[1].astype(np.int64)
    pair_dist = dist_mat[iu].astype(np.float64)

    s1, s2, tau, distance, chi, n_used = _calc_chi_spatemp_global_u_numba(
        data_raw, data_unif, pair_i, pair_j, pair_dist, tau_max, q, remove_zeros
    )

    out = pd.DataFrame({
        "s1_idx": s1,
        "s2_idx": s2,
        "s1": [colnames[i] for i in s1],
        "s2": [colnames[j] for j in s2],
        "tau": tau,
        "distance": distance,
        "chi": chi,
        "n_used": n_used
    })

    return out


#%%
###############################################################################
# 5. AGGREGATION BY DISTANCE CLASS
###############################################################################

def summarise_chi_by_distance(df_chi, dist_mat, n_bins=10):
    radius = make_radius_quantile(dist_mat, n_bins=n_bins)
    midpoints = get_midpoints(radius)

    df = df_chi.copy()
    df["dist_class"] = pd.cut(
        df["distance"],
        bins=radius,
        include_lowest=True,
        duplicates="drop"
    )

    chi_summary = (
        df.groupby(["dist_class", "tau"], observed=False)
          .agg(
              mean_chi=("chi", "mean"),
              median_chi=("chi", "median"),
              sd_chi=("chi", "std"),
              n_pairs=("chi", lambda x: np.sum(~np.isnan(x))),
              mean_n_used=("n_used", "mean"),
              min_n_used=("n_used", "min"),
              max_n_used=("n_used", "max")
          )
          .reset_index()
    )

    codes = chi_summary["dist_class"].cat.codes.to_numpy()
    valid = codes >= 0
    h = np.full(len(codes), np.nan, dtype=float)
    h[valid] = midpoints[codes[valid]]
    chi_summary["h"] = h

    return chi_summary, radius


#%%
###############################################################################
# 6. PLOTS
###############################################################################

def plot_chi_summary(chi_summary):
    plt.figure(figsize=(8, 6))

    for tau, sub in chi_summary.groupby("tau"):
        sub = sub.sort_values("h")
        plt.plot(sub["h"], sub["mean_chi"], marker="o", label=f"{tau}")

    plt.xlabel("Distance (m)")
    plt.ylabel(r"$\hat{\chi}(h,\tau)$")
    plt.legend(title="Temporal lag (5 min)", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.show()


def plot_vario_summary(chi_summary):
    df_vario = chi_summary.copy()
    df_vario["vario"] = vario_spatemp(df_vario["mean_chi"].to_numpy())

    plt.figure(figsize=(8, 6))

    for tau, sub in df_vario.groupby("tau"):
        sub = sub.sort_values("h")
        plt.plot(sub["h"], sub["vario"], marker="o", label=f"{tau}")

    plt.xlabel("Distance (m)")
    plt.ylabel(r"$\hat{\gamma}(h,\tau)$")
    plt.legend(title="Temporal lag (5 min)", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.show()

    return df_vario


def plot_chi_tau0_raw(df_chi):
    sub = df_chi[df_chi["tau"] == 0].copy()

    plt.figure(figsize=(7, 5))
    plt.scatter(sub["distance"], sub["chi"], alpha=0.5)
    plt.xlabel("Distance (m)")
    plt.ylabel(r"$\hat{\chi}(h,0)$")
    plt.tight_layout()
    plt.show()


#%%
###############################################################################
# 7. FILE PATHS
###############################################################################

filename_loc = "../../../phd_extremes/data/omsev/pluvio_mtp_loc_till_2022.csv"
filename_rain = "../../../phd_extremes/data/omsev/omsev_5min/rain_mtp_5min_2019_2025.csv"


#%%
###############################################################################
# 8. LOCATION
###############################################################################

location_gauges = pd.read_csv(filename_loc)

# adapte ici si besoin
location_gauges["Station"] = [
    "iem", "mse", "poly", "um", "cefe", "cnrs",
    "crbm", "archiw", "archie", "um35", "chu1",
    "chu2", "chu3", "chu4", "chu5", "chu6", "chu7"
]

dist_mat_full = get_dist_mat_from_lonlat(
    location_gauges,
    lon_col="Longitude",
    lat_col="Latitude"
)

print(dist_mat_full.shape)


#%%
###############################################################################
# 9. DATA
###############################################################################

rain = pd.read_csv(filename_rain)
rain_new = rain.drop(columns=["dates"], errors="ignore").copy()

for c in rain_new.columns:
    rain_new[c] = pd.to_numeric(rain_new[c], errors="coerce")

print(rain_new.shape)
print(rain_new.columns.tolist())


#%%
###############################################################################
# 10. OPTIONAL: REMOVE STATIONS IF NEEDED
###############################################################################

stations_to_remove = ["brives", "cines", "hydro"]

# enlève seulement celles qui existent
stations_to_remove = [s for s in stations_to_remove if s in rain_new.columns or s in location_gauges["Station"].values]
print("Stations to remove:", stations_to_remove)

if len(stations_to_remove) > 0:
    rain_use = rain_new.drop(columns=stations_to_remove, errors="ignore").copy()
else:
    rain_use = rain_new.copy()

# reorder location_gauges pour correspondre à l'ordre des colonnes de rain_use
location_gauges = location_gauges[~location_gauges["Station"].isin(stations_to_remove)].copy()
location_gauges = location_gauges.set_index("Station").loc[rain_use.columns].reset_index()

#%%
dist_mat = get_dist_mat_from_lonlat(location_gauges, lon_col="Longitude", lat_col="Latitude")
print("Distance matrix shape:", dist_mat.shape)

#%%
###############################################################################
# 11. CHI SPATIO-TEMPOREL
###############################################################################
#%%
from scipy.stats import rankdata
import numpy as np

site = rain_use.columns[0]
# remove nan
x_no_zeros = rain_use[site].replace(0, np.nan).to_numpy()
x_no_zeros = x_no_zeros[~np.isnan(x_no_zeros)]
x_with_zeros = rain_new[site].to_numpy()
x_with_zeros = x_with_zeros[~np.isnan(x_with_zeros)]

# -----------------------------
# 1. seuil du quantile sans zéros
q_cond_nozeros = 0.95
threshold_nozeros = np.quantile(x_no_zeros, q_cond_nozeros)

q_cond_zeros = 0.9995
threshold_zeros = np.quantile(x_with_zeros, q_cond_zeros)

print(f"Seuil quantile {q_cond_nozeros} sans zéros: {threshold_nozeros:.2f}")
print(f"Seuil quantile {q_cond_zeros} with zéros: {threshold_zeros:.2f}")

#%%

q = 0.9995
tau_max = 10
remove_zeros = False

df_chi = calc_chi_spatemp_fast_global_u(
    rain_use,
    dist_mat,
    tau_max=tau_max,
    q=q,
    remove_zeros=remove_zeros
)

print(df_chi.head())
print(df_chi.shape)
print(df_chi["chi"].describe())


#%%
###############################################################################
# 12. QUICK CHECKS
###############################################################################

print(df_chi.groupby("tau")["n_used"].describe())
plot_chi_tau0_raw(df_chi)


#%%
###############################################################################
# 13. AGGREGATION BY DISTANCE CLASS
###############################################################################

n_bins =10
chi_summary, radius = summarise_chi_by_distance(df_chi, dist_mat, n_bins=n_bins)

print("Radius breaks:", radius)
print(chi_summary.head(20))

#%%
#%%
for h in chi_summary["h"].unique()[:5]:
    sub = chi_summary[chi_summary["h"] == h]
    sub = sub.sort_values("tau")

    plt.plot(sub["tau"], sub["mean_chi"], marker="o", label=f"h={round(h)}")

plt.legend()
plt.xlabel("tau")
plt.ylabel("chi")
plt.show()

#%%
###############################################################################
# 14. PLOTS
###############################################################################

plot_chi_summary(chi_summary)
df_vario = plot_vario_summary(chi_summary)

# plot with loess smoothing
#%%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.nonparametric.smoothers_lowess import lowess

sns.set_style("whitegrid")


plt.figure(figsize=(8, 6))
for tau, sub in df_vario.groupby("tau"):
    sub = sub.sort_values("h").dropna(subset=["h", "vario"])

    # raw points
    plt.scatter(sub["h"], sub["vario"], s=5, alpha=0.7, label=f"{tau}")

    # lowess smooth
    sm = lowess(
        endog=sub["vario"].to_numpy(),
        exog=sub["h"].to_numpy(),
        frac=0.7,
        it=1,
        return_sorted=True
    )

    plt.plot(sm[:, 0], sm[:, 1], linewidth=2)

plt.xlabel("Distance (m)")
plt.ylabel(r"$\hat{\gamma}(h,\tau)$")
plt.legend(title="Temporal lag (5 min)", bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
plt.show()

#%%
###############################################################################
# 15. OPTIONAL SAVE
###############################################################################

# ex:
# out_png = "../../../phd_extremes/im/variogram/empirical/omsev/chi_spatemp_q98.png"
# plt.figure(figsize=(8, 6))
# for tau, sub in chi_summary.groupby("tau"):
#     sub = sub.sort_values("h")
#     plt.plot(sub["h"], sub["mean_chi"], marker="o", label=f"{tau}")
# plt.xlabel("Distance (m)")
# plt.ylabel(r"$\hat{\chi}(h,\tau)$")
# plt.legend(title="Temporal lag (5 min)", bbox_to_anchor=(1.02, 1), loc="upper left")
# plt.tight_layout()
# plt.savefig(out_png, dpi=300, bbox_inches="tight")
# plt.close()

#%%
#%%
import matplotlib.pyplot as plt
import seaborn as sns

df_sep = df_vario.copy()
h0 = np.sort(df_sep["h"].dropna().unique())[0]

# valeur à h0 pour chaque tau
ref = (
    df_sep.loc[df_sep["h"] == h0, ["tau", "vario"]]
    .rename(columns={"vario": "vario_h0"})
)

df_sep = df_sep.merge(ref, on="tau", how="left")
df_sep["vario_centered"] = df_sep["vario"] - df_sep["vario_h0"]

plt.figure(figsize=(8, 6))
for tau, sub in df_sep.groupby("tau"):
    sub = sub.sort_values("h")
    plt.plot(sub["h"], sub["vario_centered"], marker="o", label=f"{tau}")

plt.xlabel("Distance (m)")
plt.ylabel(r"$\gamma(h,\tau)-\gamma(h_0,\tau)$")
plt.legend(title="Temporal lag (5 min)", bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
plt.show()
# %%
#%%
tau_ref = 0
tau_cmp = 5

sub0 = df_vario[df_vario["tau"] == tau_ref][["h", "vario"]].rename(columns={"vario": "vario0"})
sub1 = df_vario[df_vario["tau"] == tau_cmp][["h", "vario"]].rename(columns={"vario": "vario1"})

df_diff = sub0.merge(sub1, on="h")
df_diff["delta"] = df_diff["vario1"] - df_diff["vario0"]

plt.figure(figsize=(7, 5))
plt.plot(df_diff["h"], df_diff["delta"], marker="o")
plt.axhline(df_diff["delta"].mean(), linestyle="--")
plt.xlabel("Distance (m)")
plt.ylabel(r"$\gamma(h,\tau_2)-\gamma(h,\tau_1)$")
plt.tight_layout()
plt.show()
# %%
