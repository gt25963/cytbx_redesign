import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Figure 2B 
# Prescreening Scores by Oligomeric State.
# C2/C3 carried forward to the full design pipeline, C4/C5 excluded.
states = ['C2', 'C3', 'C4', 'C5']
scores = [0.7616, 0.6675, 0.6305, 0.6147]
purple = '#7F77DD'
gray = '#B4B2A9'
colors = [purple, purple, gray, gray] ## first two carried forward, last two excluded

fig, ax = plt.subplots(figsize=(6, 4.5))
bars = ax.bar(states, scores, color=colors, width=0.6)
for bar, score in zip(bars, scores):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, height + 0.015,
            f'{score:.4f}', ha='center', va='bottom', fontsize=10)

ax.set_title('Prescreening Scores by Oligomeric State', fontsize=13, fontweight='bold', pad=14)
ax.set_xlabel('Oligomeric State', fontsize=11)
ax.set_ylabel('Combined Score (Boltz-2 confidence + ESM3 pTM)', fontsize=11)
ax.set_ylim(0, 0.85)

carried = mpatches.Patch(color=purple, label='Carried Forward')
excluded = mpatches.Patch(color=gray, label='Excluded (lattice incompatible)')
ax.legend(handles=[carried, excluded], loc='upper right', frameon=False, fontsize=9)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig('fig2B_prescreening.png')
print('saved to fig2B_prescreening.png')
