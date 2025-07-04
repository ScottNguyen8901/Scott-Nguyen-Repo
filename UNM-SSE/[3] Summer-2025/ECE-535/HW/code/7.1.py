import matplotlib.pyplot as plt

# Parameters
total_bw = 500  # MHz
transponder_bw = 36  # MHz
guard_bw = 4  # MHz
num_tps = 12  # from calculation

scale = 0.1  # 1 MHz = 0.1 cm (arbitrary units for matplotlib figure size)

# Calculate total used bandwidth
used_bw = num_tps * transponder_bw + (num_tps - 1) * guard_bw

fig, ax = plt.subplots(figsize=(12, 2))

# Draw total bandwidth line
ax.plot([0, total_bw*scale], [0.5, 0.5], color='black', lw=3, label='Total Bandwidth (500 MHz)')

# Draw transponders and guardbands
current_pos = 0
for i in range(num_tps):
    # Draw transponder
    ax.add_patch(plt.Rectangle((current_pos, 0.25), transponder_bw*scale, 0.5, 
                               edgecolor='blue', facecolor='skyblue', label='Transponder' if i==0 else None))
    current_pos += transponder_bw*scale
    # Draw guardband except after last transponder
    if i < num_tps - 1:
        ax.add_patch(plt.Rectangle((current_pos, 0.25), guard_bw*scale, 0.5, 
                                   edgecolor='red', facecolor='mistyrose', label='Guardband' if i==0 else None))
        current_pos += guard_bw*scale

# Label the transponders
for i in range(num_tps):
    x_center = (i * (transponder_bw + guard_bw) + transponder_bw/2) * scale
    ax.text(x_center, 0.55, f'TP{i+1}', ha='center', va='bottom', fontsize=8, color='blue')

# Format plot
ax.set_xlim(0, total_bw*scale)
ax.set_ylim(0, 1)
ax.axis('off')
ax.legend(loc='upper right')

plt.title('C-band Satellite Channeling Scheme')
plt.show()
