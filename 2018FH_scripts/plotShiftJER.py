import matplotlib.pyplot as plt
import csv

# Read data from the CSV file
data = {"mass": [], "sigma_nominal": [], "s_n_error": [], "JER_sigma_up": [], "s_up_error": [], "JER_sigma_down": [], "s_d_error": []}

with open("JER_effect.csv") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        data["mass"].append(float(row["mass"]))
        data["sigma_nominal"].append(float(row["sigma_nominal"]))
        data["s_n_error"].append(float(row["s_n_error"]))
        data["JER_sigma_up"].append(float(row["JER_sigma_up"]))
        data["s_up_error"].append(float(row["s_up_error"]))
        data["JER_sigma_down"].append(float(row["JER_sigma_down"]))
        data["s_d_error"].append(float(row["s_d_error"]))

# Create the plot
plt.figure(figsize=(10, 6))
plt.errorbar(data["mass"], data["sigma_nominal"], yerr=data["s_n_error"], fmt="o", color="black", label="sigma_nominal")
plt.errorbar(data["mass"], data["JER_sigma_up"], yerr=data["s_up_error"], fmt="o", color="blue", label="JER_sigma_up")
plt.errorbar(data["mass"], data["JER_sigma_down"], yerr=data["s_d_error"], fmt="o", color="red", label="JER_sigma_down")

# Add labels and legend
plt.xlabel("mass")
plt.ylabel("Value")
plt.title("JER Sigma Shift")
plt.legend()

# Save the plot as PNG and PDF
plt.savefig("JER_sigma_shift.png")
plt.savefig("JER_sigma_shift.pdf")

# Show the plot (optional)
plt.show()
