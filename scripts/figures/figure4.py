import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Figure 4 - panels A (C2) and B (C3)
# Chai corrected vs. aggregate scores 
c2 = {
    'cycles':     [1, 2, 3, 4],
    'seeds':      ['id96', 'id66', 'id177', 'id122'],
    'Chai real':  [0.577, 0.275, 0.187, 0.129],
    'Chai agg':   [0.688, 0.410, 0.452, 0.319],
}

c3 = {
    'cycles':     [1, 2, 3, 4, 5],
    'seeds':      ['id14', 'id35', 'id40', 'id121', 'id4'],
    'Chai real':  [0.338, 0.452, 0.580, 0.569, 0.298],
    'Chai agg':   [0.551, 0.661, 0.716, 0.718, 0.536],
}


def plot_chai_bars(data, title, outfile):
    cycles = data['cycles']
    fig, ax = plt.subplots(figsize=(7, 5))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    x = np.arange(len(cycles))
    width = 0.35

    bars1 = ax.bar(x - width / 2, data['Chai real'], width,
                    label='Chai-1 corrected (protein-protein ipTM)',
                    color='#D9822B', alpha=0.9)
    bars2 = ax.bar(x + width / 2, data['Chai agg'], width,
                    label='Chai-1 aggregate (overall score)',
                    color='#F5C089', edgecolor='#D9822B',
                    linewidth=1.2, alpha=0.9)

    for bar in list(bars1) + list(bars2):
        ax.text(bar.get_x() + bar.get_width() / 2,
                 bar.get_height() + 0.012,
                 f'{bar.get_height():.3f}',
                 ha='center', va='bottom', fontsize=7.5)

    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.set_xlabel('Cycle', fontsize=11)
    ax.set_ylabel('Chai-1 Score (ipTM)', fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels([f'Cycle {c}\n({s})'
                         for c, s in zip(cycles, data['seeds'])], fontsize=8)
    ax.set_ylim(0, 0.85)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.legend(loc='upper right', frameon=False, fontsize=9)

    plt.tight_layout()
    plt.savefig(outfile, bbox_inches='tight', dpi=300, facecolor='white', pil_kwargs={'quality': 95})
    plt.close()
    print(f"Saved {outfile}")


plot_chai_bars(c2, 'C2 Dimer: Chai-1 Corrected vs Aggregate Score', 'c2_chai_bars.jpg')
plot_chai_bars(c3, 'C3 Trimer: Chai-1 Corrected vs Aggregate Score', 'c3_chai_bars.jpg')
