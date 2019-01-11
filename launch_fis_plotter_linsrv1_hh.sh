#!/bin/sh

# fetch latest data from ingest server to Tore's home
rsync -rtvu /hs/datex/ingest/mooring/fse2/ /csys/ocean2/hhellmer/disk1/user1/work/projects/FISP/Mooring/SBD/fse2/
rsync -rtvu /hs/datex/ingest/mooring/fsw1/ /csys/ocean2/hhellmer/disk1/user1/work/projects/FISP/Mooring/SBD/fsw1/
rsync -rtvu /hs/datex/ingest/mooring/fne1/ /csys/ocean2/hhellmer/disk1/user1/work/projects/FISP/Mooring/SBD/fne1/
rsync -rtvu /hs/datex/ingest/mooring/fne2/ /csys/ocean2/hhellmer/disk1/user1/work/projects/FISP/Mooring/SBD/fne2/
wait
echo data copy complete

# Use for compiled matlab jobs, e.g. if no liscense available
# /home/csys/thatterm/fis_nrt/analysis_filchner_nrtplot_all_sites/for_redistribution_files_only/run_analysis_filchner_nrtplot_all_sites.sh /csys/ocean2/thatterm/MCRInstaller9.4

# further, none-compiled matlab tasks
#ssh -x linsrv1 "cd Dropbox/AWIsync && /opt/matlab2016a/bin/matlab -r 'mlab_job_hh'"
ssh -x linsrv1 "cd /csys/ocean2/hhellmer/disk1/user1/work/projects/FISP/Mooring/fis_nrt && /opt/matlab2016a/bin/matlab -r 'mlab_job_hh'"
wait


touch /csys/ocean2/hhellmer/disk1/user1/work/projects/FISP/Mooring/data/*.png
echo data plotted for all stations

# copy image files to infoboard
cp /csys/ocean2/hhellmer/disk1/user1/work/projects/FISP/Mooring/data/fse2_timeseries.png /csys/mob1/web-daten/oceanography/

# copy image files to webserver (could add copy of datafiles)
#scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fse2_timeseries.png thatterm@webapp-srv1a.awi.de:/var/www/apps3/FISP/
#scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fsw1_timeseries.png thatterm@webapp-srv1a.awi.de:/var/www/apps3/FISP/
#scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fne1_timeseries.png thatterm@webapp-srv1a.awi.de:/var/www/apps3/FISP/
#scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fne2_timeseries.png thatterm@webapp-srv1a.awi.de:/var/www/apps3/FISP/

echo all tasks done, bye!
exit
