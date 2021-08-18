
import numpy as np
from pygadgetreader import *
import argparse

class ReadGadget(object):
    """
    Class for reading and manipulating gadget2 files (lightcones/snapshots).

    Parameters
    ----------
    action : What to do with the files
        options: -every_N
    Notes
    -----
    returns a class for reading gadget2 files (multi or single)
    """
    def __init__(self,
                 action = "every_N"
                 ):

        self.action = action
        return

    def read_params(self):
        """
        Function to read the file parameters (lightcones/snapshots).
        Parameters
        ----------

        Return
        -----
        important parameters to read series of gadget2 files
        """
        parser = argparse.ArgumentParser(description='Reading gadget snapshots/lightcones and decrease its size')
        parser.add_argument('-kind', dest = 'kind', type=str, help='kind: snapshot or lightcone')
        parser.add_argument('-loc', dest = 'loc', type=str, help='location of the files')
        parser.add_argument('-name', dest = 'file_name', type=str, help='name: snapshot or lightcone name base - e.g., lcdm_lightcone_0, gadget2_cdm [for the lightcone it includes the number of lightcone]')
        parser.add_argument('-ini', dest = 'num_ini', type=int, help='The first file number to be used')
        parser.add_argument('-fin', dest = 'num_fin', type=int, help='The last file number to be used')
        parser.add_argument('-output', dest = 'output_name', type=str, help='output name and the address e.g., ./DATA.txt')
        assert self.action=="every_N", "The action is not chosen correctly, choose among the available options"
        if (self.action=="every_N"):
            parser.add_argument('-e', dest = 'param', type=int, help='parameter to be provided: here save every N particles only')
        args = parser.parse_args()
        # Errors
        assert args.kind and isinstance(args.kind, str) and (args.kind=="snapshot" or args.kind=="lightcone") , ' \n \033[1;31;1m Kind of the gadget2 file is not defined choose between snapshot and lightcone options e.g., \n -kind snapshot \n  for the full info run with -h option \033[0;0m'
        assert args.loc and isinstance(args.loc, str), ' \n \033[1;31;1m location of the gadget2 files is not defined use e.g., \n -loc ./DIRECTORY \n  for the full info run with -h option \033[0;0m'
        assert args.file_name and isinstance(args.file_name, str), ' \n \033[1;31;1m name of the gadget2 files is not defined use e.g., \n -name gadget2_cdm. / lcdm_lightcone \n  (use the file basis name with "." for snapshots like gadget2_cdm.xxx and for the lightcones NAMExxx_cdm "_cdm" at the end is added automatically) for the full info run with -h option \033[0;0m'
        assert isinstance(args.num_ini, int) and args.num_ini>=-1, ' \n \033[1;31;1m An integer number for the first gadget2 file is not provided use e.g., \n -ini 10 / 12 \n  (If there is only one snapshot use -i -1)  for the full info run with -h  \033[0;0m,'
        assert args.num_fin and isinstance(args.num_fin, int) and args.num_fin>=-1, ' \n \033[1;31;1m An integer number for the last gadget2 file is not provided use e.g., \n -fin 10 / 012 \n  (If there is only one snapshot use -f -1)  for the full info run with -h  \033[0;0m,'
        assert args.output_name and isinstance(args.output_name, str), ' \n \033[1;31;1m The name and format (txt/dat) of the final output  file is not defined use e.g., \n -output ./Result.txt \n  for the full info run with -h option \033[0;0m'
        if (self.action=="every_N"):
            assert args.param and isinstance(args.param, int), ' \n \033[1;31;1m An integer to decide every how many pcls to be chosen  e.g., \n -e 1000 / 012 \n  for the full info run with -h  \033[0;0m,'

        return (args.kind, args.loc, args.file_name, args.num_ini, args.num_fin, args.output_name, args.param)


    def choose_every_N(self, kind, file_location, file_base, num_ini, num_fin, output_name, save_every):
        """
        Function to read the file parameters (lightcones/snapshots) and choose every N pcls information
        ----------

        Return
        -----
        It produces the text file contaning the pcls information
        """
        print(" \n \033[1;44;30m Reading files from " + file_location+ "\033[0;0m")
        f = open (output_name,"w")
        f.write("# ID, posx, posy, posz, vx, vy, vz \n \n")
        if (kind=="lightcone"):
            print("lightcone data in gadget format requested")
            if (num_ini<0 or num_fin<0):
                print(" \n \033[1;31;1m You only chose 1 lightcone option by choosing a negative value for -ini or -fin \033[0;0m,")
                data_pos = readsnap(file_location+file_base+"_cdm","pos","dm")
                data_ID = readsnap(file_location+file_base+"_cdm","pid","dm")
                data_vel = readsnap(file_location+file_base+"_cdm","vel","dm")
                data_pos =data_pos[:,:] /1000.; # writing in Mpc/h

                for j in range (0, np.shape(data_ID)[0], save_every):
                    f.write("%d %10f %10f %10f %10f %10f %10f \r\n" %(data_ID[j],data_pos[j,0],data_pos[j,1],data_pos[j,2],
                                                                          data_vel[j,0],data_vel[j,1],data_vel[j,2]))
                    print("file number",1,"is done!")
                f.close()

            elif (num_ini>=0 and num_fin>=0):
                for i in range(num_ini,num_fin+1):

                    data_pos = readsnap(file_location+file_base+str(i)+"_cdm","pos","dm")
                    data_ID = readsnap(file_location+file_base+str(i)+"_cdm","pid","dm")
                    data_vel = readsnap(file_location+file_base+str(i)+"_cdm","vel","dm")
                    data_pos =data_pos[:,:] /1000.; # writing in Mpc/h

                    for j in range (0, np.shape(data_ID)[0], save_every):
                        f.write("%d %10f %10f %10f %10f %10f %10f \r\n" %(data_ID[j],data_pos[j,0],data_pos[j,1],data_pos[j,2],
                                                                          data_vel[j,0],data_vel[j,1],data_vel[j,2]))
                    print("file number",i,"is done!")
                f.close()


        elif (kind =="snapshot"):
            print("snapshot data in gadget format requested")
            if (num_ini<0 or num_fin<0):
                print(" \n \033[1;31;1m You only chose 1 snapshot option by choosing a negative value for -i or -f \033[0;0m,")
                data_pos = readsnap(file_location+file_base,"pos","dm")
                data_ID = readsnap(file_location+file_base,"pid","dm")
                data_vel = readsnap(file_location+file_base,"vel","dm")
                data_pos =data_pos[:,:] /1000.; # writing in Mpc/h
                for j in range (0, np.shape(data_ID)[0], save_every):
                    f.write("%d %10f %10f %10f %10f %10f %10f \r\n" %(data_ID[j],data_pos[j,0],data_pos[j,1],data_pos[j,2],
                                                                          data_vel[j,0],data_vel[j,1],data_vel[j,2]))
                print("file number",0,"is done!")
                f.close()

            elif (num_ini>=0 and num_fin>=0):
                for i in range(num_ini,num_fin+1):
                    if(i==0):
                        data_pos = readsnap(file_location+file_base+str(i),"pos","dm")[:readheader(file_location+file_base+str(1),'header')['npartThisFile'][1]]
                        data_ID = readsnap(file_location+file_base+str(i),"pid","dm")[:readheader(file_location+file_base+str(1),'header')['npartThisFile'][1]]
                        data_vel = readsnap(file_location+file_base+str(i),"vel","dm")[:readheader(file_location+file_base+str(1),'header')['npartThisFile'][1]]
                    else:
                        data_pos = readsnap(file_location+file_base+str(i),"pos","dm")[:readheader(file_location+file_base+str(i),'header')['npartThisFile'][1]]
                        data_ID = readsnap(file_location+file_base+str(i),"pid","dm")[:readheader(file_location+file_base+str(i),'header')['npartThisFile'][1]]
                        data_vel = readsnap(file_location+file_base+str(i),"vel","dm")[:readheader(file_location+file_base+str(i),'header')['npartThisFile'][1]]
                    data_pos =data_pos[:,:] /1000.; # writing in Mpc/h

                    for j in range (0, np.shape(data_ID)[0], save_every):
                        f.write("%d %10f %10f %10f %10f %10f %10f \r\n" %(data_ID[j],data_pos[j,0],data_pos[j,1],data_pos[j,2],
                                                                          data_vel[j,0],data_vel[j,1],data_vel[j,2]))
                    print("file number",i,"is done!")
                f.close()
        return

# read_file = ReadGadget("every_N");
# kind, file_location, file_base, num_ini, num_fin, output_name, param = read_file.read_params()
# read_file.choose_every_N(kind, file_location, file_base, num_ini, num_fin, output_name, param);

# Run examples:
# python read_gadget.py -kind snapshot -loc ./../output/gadget_snap/ -name gadget2_cdm -ini 0 -fin -1 -output ./txt.fine.txt -e 1
# python read_gadget.py -kind snapshot -loc ./../output/gadget_snap/ -name gadget2_cdm. -ini 0 -fin 3 -output ./txt.fine.txt -e 1
# python read_gadget.py -kind lightcone -loc ./../output/lightcones/ -name lcdm_lightcone_0  -ini 114 -fin 116 -output ./lightcone.fine.txt -e 1
