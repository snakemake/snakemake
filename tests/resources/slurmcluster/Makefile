# taken and modified from https://github.com/stevekm/slurm-cluster-vagrant

SHELL:=/bin/bash

# Create Vagrant VMs
# copy munge authentication key from slurmmaster to node
# !! need cp -p or else munge keys do not work
setup:
	vagrant up #&& \
	#vagrant up node-1 && \
	#vagrant up node-2 
	#vagrant ssh slurmmaster -- -t 'ssh-keygen -b 2048 -t rsa -q -N "" -f /home/vagrant/.ssh/id_rsa' && \
	#vagrant ssh slurmmaster -- -t 'cp /home/vagrant/.ssh/id_rsa.pub /vagrant/id_rsa.slurmmaster.pub' && \
	#vagrant ssh node-1 -- -t 'ssh-keygen -b 2048 -t rsa -q -N "" -f /home/vagrant/.ssh/id_rsa' && \
	#vagrant ssh node-1 -- -t 'cp /home/vagrant/.ssh/id_rsa.pub /vagrant/id_rsa.server.pub' && \
	#vagrant ssh node-2 -- -t 'ssh-keygen -b 2048 -t rsa -q -N "" -f /home/vagrant/.ssh/id_rsa' && \
	#vagrant ssh node-2 -- -t 'cp /home/vagrant/.ssh/id_rsa.pub /vagrant/id_rsa.server.pub' && \
	#vagrant ssh slurmmaster -- -t 'cat /vagrant/id_rsa.node-1.pub >> .ssh/authorized_keys' && \
	#vagrant ssh slurmmaster -- -t 'cat /vagrant/id_rsa.node-2.pub >> .ssh/authorized_keys' && \
	#rm -f  id_rsa.slurmmaster.pub id_rsa.server.pub

# make sure 'slurm' dir is writable for VMs
# start munge in both VMs
# start slurmctld, wait many seconds for it to fully start
# start slurmd
#start:
#	find slurm -type d -exec chmod a+rwx {} \; && \
#	vagrant ssh slurmmaster -- -t 'sudo /etc/init.d/munge start; sleep 5' && \
#	vagrant ssh node-1 -- -t 'sudo /etc/init.d/munge start; sleep 5' && \
#	vagrant ssh node-2 -- -t 'sudo /etc/init.d/munge start; sleep 5' && \
#	vagrant ssh slurmmaster -- -t 'sudo slurmctld; sleep 30' && \
#	vagrant ssh node-1 -- -t 'sudo slurmd; sleep 30' && \
#	vagrant ssh node-2 -- -t 'sudo slurmd; sleep 30' && \
#	vagrant ssh slurmmaster -- -t 'sudo scontrol update nodename=server state=resume; sinfo; sleep 5'
#
sinfo:
	vagrant ssh slurmmaster -- -t 'sinfo'

# might need this to fix node down state?
# fix:
# 	vagrant ssh slurmmaster -- -t 'sudo scontrol update nodename=server state=resume'

# https://slurm.schedmd.com/troubleshoot.html
# munge log: /var/log/munge/munged.log
test:
	@printf ">>> Checking munge keys on both machines\n"
	@vagrant ssh slurmmaster -- -t 'sudo md5sum /etc/munge/munge.key; ls -l /etc/munge/munge.key'
	@vagrant ssh server -- -t 'sudo md5sum /etc/munge/munge.key; ls -l /etc/munge/munge.key'
	@printf "\n\n>>> Checking if slurmmaster can contact node (network)\n"
	@vagrant ssh slurmmaster -- -t 'ping 10.10.10.4 -c1'
	@printf "\n\n>>> Checking if SLURM slurmmaster is running\n"
	@vagrant ssh slurmmaster -- -t 'scontrol ping'
	@printf "\n\n>>> Checking if slurmctld is running on slurmmaster\n"
	@vagrant ssh slurmmaster -- -t 'ps -el | grep slurmctld'
	@printf "\n\n>>> Checking cluster status\n"
	@vagrant ssh slurmmaster -- -t 'sinfo'
	@printf "\n\n>>> Checking if node can contact slurmmaster (network)\n"
	@vagrant ssh node-1 -- -t 'ping 10.10.10.3 -c1'
	@printf "\n\n>>> Checking if node can contact SLURM slurmmaster\n"
	@vagrant ssh node-1 -- -t 'scontrol ping'
	@printf "\n\n>>> Checking if slurmd is running on node\n"
	@vagrant ssh node-1 -- -t 'ps -el | grep slurmd'
	@printf "\n\n>>> Running a test job\n"
	@vagrant ssh slurmmaster -- -t 'sbatch --wrap="hostname"'
	@printf "\n\n>>> Running another test job\n"
	@vagrant ssh slurmmaster -- -t 'sbatch /vagrant/job.sh'
	@printf "\n\n>>> Checking node status\n"
	@vagrant ssh slurmmaster -- -t 'scontrol show nodes=server'

# pull the plug on the VMs
stop:
	vagrant halt --force slurmmaster
	vagrant halt --force node-1
	vagrant halt --force node-2

# delete the VMs
remove:
	vagrant destroy slurmmaster
	vagrant destroy node-1
	vagrant destroy node-2

# location of the SLURM default config generators for making new conf files
get-config-html:
	vagrant ssh slurmmaster -- -t 'cp /usr/share/doc/slurmctld/*.html /vagrant/'

# get rid of the SLURM log files
clean:
	find slurm -type f ! -name ".gitkeep" -exec rm -f {} \;
