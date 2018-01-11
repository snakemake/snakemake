.PHONY: run build touch stop test

build:
	docker build -t irods-server .

run: build
	docker run \
		-d \
		-p 1247:1247 \
		--name provider \
		irods-server \
			-i run_irods
	sleep 10
	docker exec -u irods provider \
		iput /incoming/infile

touch:
	docker exec -u irods provider \
		iput /incoming/infile -f

stop:
	docker ps -a \
	| grep provider \
	| cut -d ' ' -f 1 \
	| xargs -I{} \
		bash -c "docker stop {}; docker rm {}"

ls:
	docker exec -u irods provider ils -l
