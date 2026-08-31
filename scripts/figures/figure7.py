import matplotlib.pyplot as plt
import numpy as np

# Figure 7 
# "Hem1 Site" = mean distance to His37/His95
# "Hem2 Site" = mean distance to His9/His67 

structures = ['Original Scaffold\n(Design Intent)', 'AF3-Predicted\nStructure', 'Boltz-Predicted\nStructure']

dist_to_hemB_site = [
    (5.87 + 3.35) / 2, ## original scaffold: His37, His95
    (16.04 + 14.89) / 2, ## AF3: His37, His95
    (16.45 + 17.63) / 2, ## Boltz: His37, His95
]

dist_to_hemC_site = [
    (15.72 + 15.43) / 2, ## original scaffold: His9, His67
    (4.29 + 7.84) / 2, ## AF3: His9, His67
    3.81, ## Boltz: His9 only (His67 not measured)
]

x = np.arange(len(structures))
width = 0.35

fig, ax = plt.subplots(figsize=(9, 6))

light_blue = "#7FB3E8"
dark_blue = "#1B4F8C"

bars1 = ax.bar(x - width/2, dist_to_hemB_site, width,
                label='Distance to Hem1 (His37/His95)', color=dark_blue,
                edgecolor='white', linewidth=0.8)
bars2 = ax.bar(x + width/2, dist_to_hemC_site, width,
                label='Distancet to Hem2 (His9/His67)', color=light_blue,
                edgecolor='white', linewidth=0.8)

# Coordination range line
ax.axhline(y=6, color='grey', linestyle='--', linewidth=1, alpha=0.6)
ax.annotate('Approx. Coordination Range', xy=(1.0, 6), xycoords=('axes fraction', 'data'),
            xytext=(10, 0), textcoords='offset points',
            fontsize=9, color='grey', ha='left', va='center', style='italic')

ax.set_ylabel('FMN Head-Group Distance (Å)', fontsize=12, fontweight='bold')
fig.suptitle('FMN Placement: Design Intent vs. Structure-Prediction Output',
              fontsize=13, fontweight='bold', x=0.5, y=0.98)
ax.set_xticks(x)
ax.set_xticklabels(structures, fontsize=10.5)

# Legend
ax.legend(fontsize=9.5, loc='upper left', bbox_to_anchor=(1.02, 1.0), frameon=False)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(0.8)
ax.spines['bottom'].set_linewidth(0.8)
ax.set_ylim(0, 19)
ax.grid(axis='y', linestyle=':', alpha=0.3)
ax.set_axisbelow(True)

for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}', xy=(bar.get_x() + bar.get_width()/2, height),
                    xytext=(0, 4), textcoords="offset points", ha='center',
                    fontsize=9, fontweight='medium')

plt.tight_layout()
plt.savefig('fig7_placement_bias_v2.png', dpi=300, bbox_inches='tight')
print("Saved to fig7_placement_bias_v2.png")
