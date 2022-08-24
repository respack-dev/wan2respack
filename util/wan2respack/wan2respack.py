import argparse
from wan2respack_pre import preprocess
from wan2respack_core import write_wan_info
import init

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RESPACK interface code for Wannier90")  # TODO: add message
    parser.add_argument('conf', type=str, metavar="conf.toml", help='configure toml file.')
    parser.add_argument("-pp", "--pp", action="store_true",
                        help="flag of preprocess")
    args = parser.parse_args()
    conf = init.variables(args)

    if args.pp is True:
        print("Start: Run wan2respack_pre.")
        preprocess(conf)
        print("Finish: Run wan2respack_pre.")
    else:
        print("Start: Run wan2respack_core.")
        write_wan_info(conf)
        print("Finish: Run wan2respack_core.")
