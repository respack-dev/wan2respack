import os


def variables(args):
    """
    Read input arguments and a toml file (if exists) and set variables
    Returns:
        dict: conf
    """
    # Reading toml file
    if os.path.exists(args.conf):
        try:
            import tomli
            print(f"Reading {args.conf}")
            with open(args.conf, "rb") as f:
                conf = tomli.load(f)
        except:
            raise RuntimeError("Failed to read toml file.")
    else:
        raise RuntimeError("Failed to read toml file.")

    _check_essential_kwds(conf, ["base"])
    _check_essential_kwds(conf["base"], ["seedname", "QE_output_dir"])

    selfk_flag = conf["base"].get("selfk", False)
    if type(selfk_flag) is not type(True):
        raise RuntimeError("selfk in base section must be True or False")
    conf["base"]["selfk"] = selfk_flag
    if selfk_flag is True:
        if "pre" in conf.keys():
            raise RuntimeError("pre section must not exist when selfk is True")
    else:
        _check_essential_kwds(conf["pre"], ["ref", "output"])
        for section_name in ["ref", "output"]:
            try:
                _check_essential_kwds(conf["pre"][section_name], ["nscf", "win"])
            except:
                raise RuntimeError("Setting of {} section.".format(section_name))

    conf["gen_mk"] = os.path.join(os.path.dirname(__file__), "gen_mk.x")
    conf["gen_wan"] = os.path.join(os.path.dirname(__file__), "gen_wan.x")
    _check_files([conf["base"]["QE_output_dir"], conf["gen_mk"], conf["gen_wan"]])

    return conf

def _check_essential_kwds(dict, check_kwds):
    for kwd_name in check_kwds:
        if kwd_name not in dict.keys():
            raise RuntimeError("{} is not found.".format(kwd_name))

def _check_files(check_files):
    for file_name in check_files:
        if not os.path.exists(file_name):
            raise RuntimeError("{} is not found.".format(file_name))
