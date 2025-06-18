import matplotlib.pyplot as plt
import statistics

with open("agg_data_B-D-Glc_OPC.txt") as f:
    contents = f.read().splitlines()
    
frames = list(range(1,1001))
agg_num = []
avg_agg_size = []

for line in contents:
    agg_num.append(float(line.split()[0]))
    avg_agg_size.append(float(line.split()[1]))

fig,ax = plt.subplots(2)

ax[0].plot(frames, agg_num, color = "green")
ax[0].set_xlabel("Frame")
ax[0].set_ylabel("Number of Aggregates")

ax[1].plot(frames, avg_agg_size, color = "purple")
ax[1].set_xlabel("Frame")
ax[1].set_ylabel("Average Aggregate Size")

fig.tight_layout()
fig.savefig("B-D-Glc_OPC_agg_over_time.png")

# Now make a bar graph with stdevs for the average number of aggregates and average
# of the average aggregate size. 


fig,ax = plt.subplots(1)

ax.bar(["Average # of Aggregates", "Average Aggregate Size"],[statistics.mean(agg_num), statistics.mean(avg_agg_size)], width = 0.5, yerr = [statistics.stdev(agg_num), statistics.stdev(avg_agg_size)], capsize = 5 , color = ["green","purple"])
ax.set_ylim(0,5)
ax.set_title("Averaged over Full Simulation")

fig.tight_layout()
fig.savefig("B-D-Glc_OPC_agg_data_bar_graph.png")

print(round(statistics.mean(agg_num), 2), round(statistics.stdev(agg_num),2))
print(round(statistics.mean(avg_agg_size),2), round(statistics.stdev(avg_agg_size),2))