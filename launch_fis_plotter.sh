#!/bin/sh

# fetch latest data from ingest server to Tore's home
rsync -rtvu /hs/datex/ingest/mooring/fse2/ /home/csys/thatterm/Dropbox/AWIsync/SBD/fse2/
rsync -rtvu /hs/datex/ingest/mooring/fsw1/ /home/csys/thatterm/Dropbox/AWIsync/SBD/fsw1/
rsync -rtvu /hs/datex/ingest/mooring/fne1/ /home/csys/thatterm/Dropbox/AWIsync/SBD/fne1/
rsync -rtvu /hs/datex/ingest/mooring/fne2/ /home/csys/thatterm/Dropbox/AWIsync/SBD/fne2/
rsync -rtvu /hs/datex/ingest/mooring/fse2/ /home/csys/thatterm/SBD/
rsync -rtvu /hs/datex/ingest/mooring/fsw1/ /home/csys/thatterm/SBD/
rsync -rtvu /hs/datex/ingest/mooring/fne1/ /home/csys/thatterm/SBD/
rsync -rtvu /hs/datex/ingest/mooring/fne2/ /home/csys/thatterm/SBD/
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

# plotting
# old command reading data from Tore's home
#/home/csys/thatterm/fis_nrt/analysis_filchner_nrtfun_plot/for_redistribution_files_only/run_analysis_filchner_nrtfun_plot.sh /home/csys/thatterm/v901 /home/csys/thatterm/Dropbox/Osci/FISP/#DATA/SBD/300234061032780_ fse2 /home/csys/thatterm/fis_nrt /csys/mob1/web-daten/oceanography/

# new plotting routines reading data in ingest archive
#/home/csys/thatterm/fis_nrt/analysis_filchner_nrtfun_plot/for_redistribution_files_only/run_analysis_filchner_nrtfun_plot.sh /home/csys/thatterm/v901 /home/csys/thatterm/SBD/300234061032780_ fse2 /home/csys/thatterm/fis_nrt /home/csys/thatterm/fis_nrt/analysis_filchner_nrtfun_plot/
#wait
#echo plotted fse2

#/home/csys/thatterm/fis_nrt/analysis_filchner_nrtfun_plot/for_redistribution_files_only/run_analysis_filchner_nrtfun_plot.sh /home/csys/thatterm/v901 /home/csys/thatterm/SBD/300234061031800_ fsw1 /home/csys/thatterm/fis_nrt /home/csys/thatterm/fis_nrt/analysis_filchner_nrtfun_plot/
#wait
#echo plotted fsw1

/home/csys/thatterm/fis_nrt/analysis_filchner_nrtplot_all_sites/for_redistribution_files_only/run_analysis_filchner_nrtplot_all_sites.sh /home/csys/thatterm/v901
wait
touch /home/csys/thatterm/Dropbox/AWIsync/SBD/*.png
echo data plotted for all stations

# copy image files to infoboard and webserver
cp /home/csys/thatterm/Dropbox/AWIsync/SBD/fse2_timeseries.png /csys/mob1/web-daten/oceanography/
scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fse2_timeseries.png thatterm@rep3-vm.awi.de:/var/www/FISP/
scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fsw1_timeseries.png thatterm@rep3-vm.awi.de:/var/www/FISP/
scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fne1_timeseries.png thatterm@rep3-vm.awi.de:/var/www/FISP/
scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fne2_timeseries.png thatterm@rep3-vm.awi.de:/var/www/FISP/

sleep 900
echo waited 900 s second time for dropbox synching

ps aux | grep -i dropbox | awk {'print $2'} | xargs kill -9
echo dropbox killed

echo all tasks done, bye!
exit
