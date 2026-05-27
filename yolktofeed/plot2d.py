import matplotlib.pyplot as plt
import numpy as np
from rcparams import rc_update
from matplotlib import rcParams

def plot_pca_components(pca_t, Xpca, sample_meta, DAYS: list):
    rc_update(rcParams)
    DAY_COLORS = {
        4:  '#2ecc71',   # green
        6:  '#27ae60',   # dark green
        8:  '#f1c40f',   # yellow
        10: '#e67e22',   # orange
        12: '#e74c3c',   # red
        14: '#8e44ad',   # PURPLE ← Day 14
        16: '#3498db',   # blue
        18: '#1abc9c',   # teal
        20: '#2c3e50',   # dark navy
        }

    fig, axes = plt.subplots(1, 1, figsize=(8,8))
    fig.suptitle('PC1 vs PC2', fontsize=14, fontweight='bold')

    ax = axes
    for day in DAYS:
        idx = np.where(sample_meta['Day'].values == day)[0]
        ax.scatter(Xpca[idx, 0], Xpca[idx, 1],c=DAY_COLORS[day], 
                        label=f'Day {day}', s=90, 
                        edgecolors='white', zorder=3)

        ax.set_xlabel(f'PC1 ({pca_t.explained_variance_ratio_[0]*100:.1f}%)', fontsize=12)
        ax.set_ylabel(f'PC2 ({pca_t.explained_variance_ratio_[1]*100:.1f}%)', fontsize=12)
        ax.set_xlim((-300,300))
        ax.set_ylim((-100,100))
        ax.set_title(f'PCA for days in {DAYS}', fontsize=12)
        ax.legend(title='Day', bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=9)
        ax.grid(True, alpha=0.4)
    
    
    plt.tight_layout()
    plt.savefig(f"results/pca_components_in{DAYS}.pdf".replace(" ",""))


def volcano_plot(de_df, up, dn):
    rc_update(rcParams)

    fig, axes = plt.subplots(1, 1, figsize=(6, 6))
    ax = axes
    ns = de_df[~de_df['sig']]
    comp_days = f'{de_df.attrs['days_b']}vs{de_df.attrs['days_a']}'
    comp_days.replace(" ","")

    ax.scatter(ns['LFC'],  ns['-log10p'],  c='#bdc3c7', s=5,  alpha=0.4, label='Not significant')
    ax.scatter(up['LFC'],  up['-log10p'],  c='#e74c3c', s=15, alpha=0.8, label=f'Up at {de_df.attrs['days_b']} ({len(up)})')
    ax.scatter(dn['LFC'],  dn['-log10p'],  c='#3498db', s=15, alpha=0.8, label=f'Down at {de_df.attrs['days_b']} ({len(dn)})')

    ax.axhline(-np.log10(0.05), color='gray', linestyle='--', linewidth=1, alpha=0.7)
    ax.axvline( 1, color='#e74c3c', linestyle='--', linewidth=1, alpha=0.4)
    ax.axvline(-1, color='#3498db', linestyle='--', linewidth=1, alpha=0.4)

    for _, row in up.nlargest(4, '-log10p').iterrows():
        ax.annotate(row['Gene'], (row['LFC'], row['-log10p']),
             fontsize=7, color='#c0392b', xytext=(3, 2), textcoords='offset points')
    for _, row in dn.nlargest(4, '-log10p').iterrows():
        ax.annotate(row['Gene'], (row['LFC'], row['-log10p']),
             fontsize=7, color='#2980b9', xytext=(3, 2), textcoords='offset points')
    ax.set_xlabel(f'Log2 Fold Change ({de_df.attrs['days_b']} vs {de_df.attrs['days_a']})', fontsize=12)
    ax.set_ylabel('-log10(p-value)', fontsize=12)
    ax.set_title('Volcano Plot', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.savefig(f"results/volcanoPlotFor{comp_days}.png")
