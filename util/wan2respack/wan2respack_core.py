import argparse
import os.path
import subprocess
from qe2respack import qe2respack
import init


def write_wan_info(conf: dict):
    """
    Generate input files for RESPACK.
    :return:
    """

    path_to_dir_wfn = "dir-wfn"
    path_to_dir_wan = "dir-wan"
    seed_name = conf['base']['seedname']
    path_to_dir_wfn_irre = "{}-irre".format(path_to_dir_wfn)
    path_to_dir_wfn_full = "{}-full".format(path_to_dir_wfn)

    for file_name in ["{}.chk".format(seed_name), "{}.win".format(seed_name), path_to_dir_wfn]:
        if not os.path.exists(file_name):
            raise RuntimeError("{} is not found.".format(file_name))

    # Rename path_to_dir_wfn by path_to_dir_wfn_irre
    # for replacing dir_wfn by original dif_wfn at the last procedure
    os.rename(path_to_dir_wfn, path_to_dir_wfn_irre)
    if os.path.exists(path_to_dir_wan):
        os.rename(path_to_dir_wan, "{}.bak".format(path_to_dir_wan))
    os.mkdir(path_to_dir_wan)

    print("Start: running qe2respack.")
    qe2respack(conf["base"]["QE_output_dir"])
    print("Finish: running qe2respack.")

    print("Start: running {}.".format(conf["gen_wan"]))
    subprocess.run([conf["gen_wan"], seed_name])
    print("Finish: running {}.".format(conf["gen_wan"]))

    # Replace dir_wfn by original dif_wfn
    os.rename(path_to_dir_wfn, path_to_dir_wfn_full)
    os.rename(path_to_dir_wfn_irre, path_to_dir_wfn)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RESPACK interface code for Wannier90")  # TODO: add message
    parser.add_argument('conf', type=str, metavar="conf.toml", help='configure toml file.')
    args = parser.parse_args()

    conf = init.variables(args)

    write_wan_info(conf)
