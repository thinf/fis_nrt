#!/bin/sh

# fetch latest data from ingest server to Tore's home
rsync -rtvu /hs/datex/ingest/mooring/fse2/ /home/csys/thatterm/Dropbox/AWIsync/SBD/fse2/
rsync -rtvu /hs/datex/ingest/mooring/fsw1/ /home/csys/thatterm/Dropbox/AWIsync/SBD/fsw1/
rsync -rtvu /hs/datex/ingest/mooring/fne1/ /home/csys/thatterm/Dropbox/AWIsync/SBD/fne1/
rsync -rtvu /hs/datex/ingest/mooring/fne2/ /home/csys/thatterm/Dropbox/AWIsync/SBD/fne2/
rsync -rtvu /hs/datex/ingest/mooring/fse2/ /home/csys/thatterm/SBD/fse2/
rsync -rtvu /hs/datex/ingest/mooring/fsw1/ /home/csys/thatterm/SBD/fsw1/
rsync -rtvu /hs/datex/ingest/mooring/fne1/ /home/csys/thatterm/SBD/fne1/
rsync -rtvu /hs/datex/ingest/mooring/fne2/ /home/csys/thatterm/SBD/fne2/
wait
echo data copy complete

## sync with dropbox to access data on Tore's local machine
echo start
~/.dropbox-dist/dropboxd &
echo dropbox started

# wait for dropbox syncing and shut down
#sleep 300
sleep 900
echo waited 900 s h for dropbox synching

#/home/csys/thatterm/fis_nrt/analysis_filchner_nrtplot_all_sites/for_redistribution_files_only/run_analysis_filchner_nrtplot_all_sites.sh /home/csys/thatterm/v901
#wait

/home/csys/thatterm/fis_nrt/analysis_filchner_nrtplot_all_sites/for_redistribution_files_only/run_analysis_filchner_nrtplot_all_sites.sh /csys/ocean2/thatterm/MCRInstaller9.4
# further, none-compiled matlab tasks
ssh -x linsrv1 "cd Dropbox/AWIsync && /opt/matlab2016a/bin/matlab -r 'mlab_job'"
wait


touch /home/csys/thatterm/Dropbox/AWIsync/SBD/*.png
echo data plotted for all stations

# copy image files to infoboard and webserver
cp /home/csys/thatterm/Dropbox/AWIsync/SBD/fse2_timeseries.png /csys/mob1/web-daten/oceanography/
#scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fse2_timeseries.png thatterm@rep3-vm.awi.de:/var/www/FISP/
#scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fsw1_timeseries.png thatterm@rep3-vm.awi.de:/var/www/FISP/
#scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fne1_timeseries.png thatterm@rep3-vm.awi.de:/var/www/FISP/
#scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fne2_timeseries.png thatterm@rep3-vm.awi.de:/var/www/FISP/

scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fse2_timeseries.png thatterm@webapp-srv1a.awi.de:/var/www/apps3/FISP/
scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fsw1_timeseries.png thatterm@webapp-srv1a.awi.de:/var/www/apps3/FISP/
scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fne1_timeseries.png thatterm@webapp-srv1a.awi.de:/var/www/apps3/FISP/
scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fne2_timeseries.png thatterm@webapp-srv1a.awi.de:/var/www/apps3/FISP/

## further, none-compiled matlab tasks
#ssh -x linsrv1 "cd Dropbox/AWIsync && /opt/matlab2016a/bin/matlab -r 'mlab_job'"
#wait

# need to find alternative password protected data storage
#scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fsw1.mat thatterm@rep3-vm.awi.de:/home/thatterm/
#scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fse2.mat thatterm@rep3-vm.awi.de:/home/thatterm/
#scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fne1.mat thatterm@rep3-vm.awi.de:/home/thatterm/
#scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fne2.mat thatterm@rep3-vm.awi.de:/home/thatterm/

sleep 900
echo waited 900 s second time for dropbox synching

ps aux | grep -i dropbox | awk {'print $2'} | xargs kill -9
echo dropbox killed

echo all tasks done, bye!
exit
