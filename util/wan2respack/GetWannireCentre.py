#!/usr/bin/env python
import argparse
import os
import numpy as np   

def main():
    parser = argparse.ArgumentParser(description="RESPACK interface code for Wannier90")  # TODO: add message
    parser.add_argument('conf', type=str, metavar="conf.toml", help='configure toml file.')
    parser.add_argument("-pp", "--pp", action="store_true",
                        help="flag of preprocess")
    args = parser.parse_args()
    conf = simple_variables(args)
    seedname = conf["base"]["seedname"]
 
    file_name = "%s.wout" % (seedname)
    vec_lat   = GetLatticeVec(file_name)
    vec_wan   = GetWanCentre(file_name)
    frac_wan  = Convert2Frac(vec_lat,vec_wan)
    #print(vec_lat)
    #print(vec_wan)
    #print(frac_wan)
    with open("zvo_geom.dat", 'w') as f:
        print("  %f %f %f " % (vec_lat[0][0],vec_lat[0][1],vec_lat[0][2]),file=f)
        print("  %f %f %f " % (vec_lat[1][0],vec_lat[1][1],vec_lat[1][2]),file=f)
        print("  %f %f %f " % (vec_lat[2][0],vec_lat[2][1],vec_lat[2][2]),file=f)
        print("  %d " % (frac_wan.shape[0]),file=f)
        for cnt in range(frac_wan.shape[0]):
            print("  %f %f %f " % (frac_wan[cnt][0],frac_wan[cnt][1],frac_wan[cnt][2]),file=f)


def Convert2Frac(vec_lat,vec_wan):
    inv_vec_lat = np.linalg.inv(vec_lat)
    num_wan     = vec_wan.shape[0]
    frac_wan     = np.zeros((num_wan,3),dtype=np.float64)
    #print(inv_vec_lat,num_wan)
    for cnt in range(num_wan):
        for idx in range(3):
            tmp = 0.0
            for cnt_i in range(3):
                tmp += vec_wan[cnt][cnt_i]*inv_vec_lat[cnt_i][idx]
            frac_wan[cnt][idx] = tmp
    return frac_wan

def GetLatticeVec(file_name):
    with open("%s" % (file_name)) as f:
        tmp      = f.read()
        tmp      = tmp.split("\n")
    for cnt in range(len(tmp)):
        tmp_2 = tmp[cnt].split()
        if len(tmp_2)> 0 and tmp_2[0] == "Lattice" and tmp_2[1]=="Vectors":
            cnt_s = cnt+1
            break
    vec_lat  = np.zeros((3,3),dtype=np.float64)
    idx      = 0
    for cnt in range(cnt_s,cnt_s+3):
        tmp_2 = tmp[cnt].split()
        vec_lat[idx][0] = float(tmp_2[1])
        vec_lat[idx][1] = float(tmp_2[2])
        vec_lat[idx][2] = float(tmp_2[3])
        idx +=1
    return vec_lat

def GetWanCentre(file_name):
    with open("%s" % (file_name)) as f:
        tmp      = f.read()
        tmp      = tmp.split("\n")
    for cnt in range(len(tmp)):
        tmp_2 = tmp[cnt].split()
        if len(tmp_2)> 0 and tmp_2[0] == "Final" and tmp_2[1]=="State":
            cnt_s = cnt+1
            #print(tmp_2[0],tmp_2[1])
            break
    cnt_e=cnt_s
    for cnt in range(cnt_s,len(tmp)):
        tmp_2 = tmp[cnt].split()
        if len(tmp_2)> 0 and tmp_2[0] == "WF" and tmp_2[1]=="centre":
            cnt_e+=1
        else:
            break
    num_wan = cnt_e-cnt_s
    vec_wan  = np.zeros((num_wan,3),dtype=np.float64)
    #print(cnt_s,cnt_e,num_wan)
    for cnt in range(cnt_s,cnt_e):
        tmp_2 = tmp[cnt].split()
        #print(tmp_2[4],tmp_2[6],tmp_2[7],tmp_2[8],tmp_2[10])
        cnt_wan = int(tmp_2[4])-1  
        vec_wan[cnt_wan][0] = float(tmp_2[6].rstrip(','))
        vec_wan[cnt_wan][1] = float(tmp_2[7].rstrip(','))
        vec_wan[cnt_wan][2] = float(tmp_2[8].rstrip(','))
    return vec_wan

def simple_variables(args):
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

    return conf
       

if __name__ == "__main__":
    main()
