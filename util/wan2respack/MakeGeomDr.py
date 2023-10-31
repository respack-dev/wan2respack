#!/usr/bin/env python
import argparse
import os
import numpy as np   
import hwave.qlms

def main():
    parser = argparse.ArgumentParser(description="RESPACK interface code for Wannier90")  # TODO: add message
    parser.add_argument('conf', type=str, metavar="conf.toml", help='configure toml file.')
    parser.add_argument("-pp", "--pp", action="store_true",
                        help="flag of preprocess")
    args = parser.parse_args()
    conf = simple_variables(args)
    seedname = conf["base"]["seedname"]
 
    #[s] make zvo_geom.dat
    file_name = "%s.wout" % (seedname)
    vec_lat   = GetLatticeVec(file_name)
    vec_wan   = GetWanCentre(file_name)
    frac_wan  = Convert2Frac(vec_lat,vec_wan)
    with open("zvo_geom.dat", 'w') as f:
        print("  %f %f %f " % (vec_lat[0][0],vec_lat[0][1],vec_lat[0][2]),file=f)
        print("  %f %f %f " % (vec_lat[1][0],vec_lat[1][1],vec_lat[1][2]),file=f)
        print("  %f %f %f " % (vec_lat[2][0],vec_lat[2][1],vec_lat[2][2]),file=f)
        print("  %d " % (frac_wan.shape[0]),file=f)
        for cnt in range(frac_wan.shape[0]):
            print("  %f %f %f " % (frac_wan[cnt][0],frac_wan[cnt][1],frac_wan[cnt][2]),file=f)
    #[e] make zvo_geom.dat
   

    #[s] make zvo_dr.dat
    name_hr = "%s_hr.dat" % (seedname)
    Lx,Ly,Lz,orb_num = read_w90(name_hr)
    print(Lx, Ly,Lz,orb_num)
    #
    Ncond = int(1*Lx*Ly*Lz*orb_num/3)
    #
    MakeInputToml("input.toml","zvo_geom.dat",name_hr,Lx,Ly,Lz,Ncond)
    #os.system("hwave input.toml")
    hwave.qlms.run(input_file="input.toml") 
 
    data  = np.load("output/green.npz")
    green = data["green"]

    name_out = "zvo_dr.dat"
    with open(name_out, "w") as fw:
        print("# zvo_dr.dat by wannier90 format",file=fw)
        print("%d " % (orb_num),file=fw)
        Nlattice = Lx*Ly*Lz
        print("%d " % (Nlattice),file=fw)
        for idx in range(Nlattice):
            print(" 1 ",end="",file=fw)
            if (idx+1)%15 == 0:
                print(" ",file=fw)
        if (idx+1)%15 != 0:
            print(" ",file=fw)
        for tmp_cnt_x in range(-int(Lx/2),int(Lx/2)+1,1):
            cnt_x = ConvertCnt(Lx,tmp_cnt_x)
            for tmp_cnt_y in range(-int(Ly/2),int(Ly/2)+1,1):
                cnt_y = ConvertCnt(Ly,tmp_cnt_y)
                for tmp_cnt_z in range(-int(Lz/2),int(Lz/2)+1,1):
                    cnt_z = ConvertCnt(Lz,tmp_cnt_z)
                    cnt_tot = cnt_z+Lz*cnt_y+Lz*Ly*cnt_x
                    for orb_j in range(orb_num):
                        for orb_i in range(orb_num):
                            tmp_val  = green[cnt_tot][0][orb_i][0][orb_j]
                            tmp_val += green[cnt_tot][1][orb_i][1][orb_j]
                            print(" %d %d %d %d %d %f %f " % (tmp_cnt_x,tmp_cnt_y,tmp_cnt_z,orb_i+1,orb_j+1,tmp_val.real,tmp_val.imag),file=fw)


def MakeInputToml(name_in,name_geom,name_hr,Lx,Ly,Lz,Ncond):
    with open(name_in, 'w') as fw:
        print("[log]",file=fw)
        print("  print_level = 1",file=fw)
        print("  print_step = 10",file=fw)
        print("[mode]",file=fw)
        print("  mode = \"UHFk\"  ",file=fw)
        print("[mode.param]",file=fw)
        print(" # 2Sz = 0",file=fw)
        print("  Ncond = %d" % (Ncond),file=fw)
        print("  IterationMax = 1000",file=fw)
        print("  EPS = 8",file=fw)
        print("  Mix = 0.5",file=fw)
        print("  T = 0.0",file=fw)
        print("  CellShape = [ %d, %d, %d ]" % (Lx,Ly,Lz),file=fw)
        print("  SubShape = [ 1, 1, 1 ]",file=fw)
        print("[file]",file=fw)
        print("[file.input.interaction]",file=fw)
        print("  path_to_input = \"./\"",file=fw)
        print("  Geometry = \"%s\"" % (name_geom),file=fw)
        print("  Transfer = \"%s\"" % (name_hr),file=fw)
        print("[file.output]",file=fw)
        print("  path_to_output = \"output\"",file=fw)
        print("  energy = \"energy.dat\"",file=fw)
        print("  eigen = \"eigen\"",file=fw)
        print("  green = \"green\"",file=fw)
    
def read_w90(name_in):
    with open(name_in, 'r') as f:
        l_strip = [s.strip() for s in f.readlines()[1:]]

    nr = int(l_strip[1])
    nints_per_line = 15
    skip_line = nr // nints_per_line
    if nr % nints_per_line != 0:
        skip_line += 1
    L_max   = np.zeros((4),dtype=np.int64)
    L_min  = np.zeros((4),dtype=np.int64)
    for idx, line in enumerate(l_strip[2 + skip_line:]):
        values = line.split()
        # if data is empty, break
        if len(values) == 0:
            break
        for cnt in range(len(L_max)):
            if int(values[cnt]) > L_max[cnt]:
                L_max[cnt] = int(values[cnt])
            if int(values[cnt]) < L_min[cnt]:
                L_min[cnt] = int(values[cnt])


    return L_max[0]-L_min[0]+1,L_max[1]-L_min[1]+1,L_max[2]-L_min[2]+1,L_max[3]-L_min[3]

def ConvertCnt(max_L,cnt):
    if cnt<0:
        cnt = cnt+max_L
    return cnt

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
