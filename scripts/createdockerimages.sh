echo -e "-----> Building trimAl Docker Image\n"
docker build -t trimal -f DockerFolder/Dockerfile.trimAl .
echo -e "\n\n-----> Building readAl Docker Image\n"
docker build -t readal -f DockerFolder/Dockerfile.readAl .
echo -e "\n\n-----> Building statAl Docker Image\n"
docker build -t statal -f DockerFolder/Dockerfile.statAl .
