import argparse
import os.path
import subprocess
import shutil
import init
from qe2respack import qe2respack

def replace(mode: str, ref: str, new: str):
    if not os.path.exists(ref):
        raise RuntimeError(f"{ref} is not found.")
    if os.path.exists(new):
        os.rename(new, new + ".bak")

    if mode == "nscf":
        replace_nscf(ref, new)
    elif mode == "win":
        replace_win(ref, new)


def replace_nscf(ref_nscf: str, new_nscf: str):
    with open(ref_nscf) as fp:
        lines = fp.readlines()

    new_lines = ""
    for line in lines:
        if "K_POINTS" in line:
            break
        new_lines += line
        if "&system" in line:
            new_lines += " nosym = .true.\n"

    new_lines += "K_POINTS {crystal}\n"
    with open("dat.sample_mk") as fp:
        new_lines += fp.readline()
        new_lines += "".join([line.strip() + "  1\n" for line in fp.readlines()])

    with open(new_nscf, "w") as fp:
        fp.write(new_lines)


def replace_win(ref_win: str, new_win: str):
    shutil.copy(ref_win, new_win)

    kp_lines = "begin kpoints\n"
    with open("dat.sample_mk") as fp:
        fp.readline()
        kp_lines += "".join(fp.readlines())
    kp_lines += "end kpoints\n"

    with open(new_win, "a") as fp:
        fp.write(kp_lines)


def preprocess(conf: dict):
    # Export "dat.sample_mk" and "dat.kg_respack"

    print("Start: Run qe2respack")
    if os.path.exists('dir-wfn'):
        print('dir-wfn/ already exists. qe2respack.py moves this to dir-wfn.backup/ as a backup.')
        if os.path.exists('dir-wfn.backup'):
            shutil.rmtree('dir-wfn.backup')
        shutil.move('dir-wfn', 'dir-wfn.backup')
    qe2respack(conf["base"]["QE_output_dir"])
    print("End: Run qe2respack")

    print("Start: Run gen_mk")
    subprocess.run([conf["gen_mk"]])
    print("End: Run gen_mk")

    for word in ["nscf", "win"]:
        if conf["base"]["selfk"] is False:
            print("Generating {} input for Wannier90, {}, from {}."
                  .format(word, conf["pre"]["output"][word], conf["pre"]["ref"][word]))
            replace(word, conf["pre"]["ref"][word], conf["pre"]["output"][word])
        else:
            print("Please use dat.sample_mk for k points in {}.".format(word))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RESPACK interface code for Wannier90")
    parser.add_argument('conf', type=str, metavar="conf.toml", help='configure toml file.')
    args = parser.parse_args()
    conf = init.variables(args)

    preprocess(conf)
