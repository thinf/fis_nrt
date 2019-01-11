#!/bin/sh

# fetch latest data from ingest server to Tore's home
rsync -rtvu /hs/datex/ingest/mooring_onisi/fse2/ /home/csys/thatterm/Dropbox/AWIsync/SBD/fse2/
rsync -rtvu /hs/datex/ingest/mooring_onisi/fsw1/ /home/csys/thatterm/Dropbox/AWIsync/SBD/fsw1/
rsync -rtvu /hs/datex/ingest/mooring_onisi/fne1/ /home/csys/thatterm/Dropbox/AWIsync/SBD/fne1/
rsync -rtvu /hs/datex/ingest/mooring_onisi/fne2/ /home/csys/thatterm/Dropbox/AWIsync/SBD/fne2/
rsync -rtvu /hs/datex/ingest/mooring_onisi/fse2/ /home/csys/thatterm/SBD/fse2/
rsync -rtvu /hs/datex/ingest/mooring_onisi/fsw1/ /home/csys/thatterm/SBD/fsw1/
rsync -rtvu /hs/datex/ingest/mooring_onisi/fne1/ /home/csys/thatterm/SBD/fne1/
rsync -rtvu /hs/datex/ingest/mooring_onisi/fne2/ /home/csys/thatterm/SBD/fne2/
wait
echo data copy complete

# further, none-compiled matlab tasks
ssh -x linsrv1 "cd Dropbox/AWIsync && /opt/matlab2016a/bin/matlab -r 'mlab_job'"
wait

touch /home/csys/thatterm/Dropbox/AWIsync/SBD/*.png
echo data plotted for all stations

# copy image files to infoboard and webserver
cp /home/csys/thatterm/Dropbox/AWIsync/SBD/fse2_timeseries.png /csys/mob1/web-daten/oceanography/

scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fse2_timeseries.png thatterm@webapp-srv1a.awi.de:/var/www/apps3/FISP/
scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fsw1_timeseries.png thatterm@webapp-srv1a.awi.de:/var/www/apps3/FISP/
scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fne1_timeseries.png thatterm@webapp-srv1a.awi.de:/var/www/apps3/FISP/
scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fne2_timeseries.png thatterm@webapp-srv1a.awi.de:/var/www/apps3/FISP/

# use Dropbox uploader to upload new files
# https://github.com/andreafabrizi/Dropbox-Uploader 
~/Dropbox-Uploader/dropbox_uploader.sh -s upload ~/Dropbox/AWIsync/* ./
## further, none-compiled matlab tasks
#ssh -x linsrv1 "cd Dropbox/AWIsync && /opt/matlab2016a/bin/matlab -r 'mlab_job'"
#wait

# need to find alternative password protected data storage
#scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fsw1.mat thatterm@rep3-vm.awi.de:/home/thatterm/
#scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fse2.mat thatterm@rep3-vm.awi.de:/home/thatterm/
#scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fne1.mat thatterm@rep3-vm.awi.de:/home/thatterm/
#scp /home/csys/thatterm/Dropbox/AWIsync/SBD/fne2.mat thatterm@rep3-vm.awi.de:/home/thatterm/


#ps aux | grep -i dropbox | awk {'print $2'} | xargs kill -9
#echo dropbox killed

echo all tasks done, bye!
exit
