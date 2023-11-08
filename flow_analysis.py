from flow_analysis.readers.read_hirep import read_flows_hirep
from flow_analysis.measurements.scales import measure_w0
import csv, os

os.makedirs("output/Qhistories/", exist_ok=True)
path = 'input_data/wilson_flow/'


def process_flow_file_hirep_for_singlets(hirep_file, W0, filename):
    # obtain topological charge Q, the number of the configurations as specified in the configuration filenames and the
    # gradient flow scale. For the ensembles in this dataset the reference scale has been fixed to W0=0.28125 in the
    # measurements.
    flows = read_flows_hirep(hirep_file)
    Qs = flows.Q_history()
    trajectories = flows.trajectories
    scale = measure_w0(flows, W0)

    # Write everything into a csv with that is compatible with the current julia scripts used in the analysis
    f = open(filename, "w")
    f.write("trajectory,Q (w0 = %s)\n" % scale)
    for i in zip(trajectories, Qs):
        f.write("%s,%s\n" % (i[0], i[1]))
    f.close()


with open("input/parameters/param_flow.csv", newline="") as prmfile:
    filereader = csv.reader(prmfile, delimiter=";")
    next(filereader, None)  # skip the headers
    for row in filereader:
        flow_file = row[0]
        hirep_file = path + flow_file + "/out_flow"
        W0 = 0.28125
        process_flow_file_hirep_for_singlets(
            hirep_file, W0, "output/Qhistories/" + os.path.basename(flow_file)
        )
