## Already-built PowerBacGWAS Docker images

You can use already built Docker images for PowerBacGWAS for MacOS arm64 and Linux amd64 platforms, available on https://hub.docker.com/repository/docker/francesccoll/powerbacgwas, under the following names: 'francesccoll/powerbacgwas:arm64' and 'francesccoll/powerbacgwas:amd64', respectively.

If you need to run the Docker image in a different host architecture, see steps below on how to build the Docker image.

## Steps to build PowerBacGWAS Docker image

* Create an empty directory
```console
mkdir docker_image
```
* Download all PowerBacGWAS scripts from GitHub and copy them into this empty directory
```console
cd docker_image
git clone https://github.com/francesccoll/powerbacgwas
cp ./powerbacgwas/scripts/* .
```
* Copy the Dockerfile into this directory and delete the downloaded GitHub directory
```console
cp ./powerbacgwas/docker/Dockerfile .
rm -fr powerbacgwas
```
* Build the Docker image for your platform of choice
```console
docker buildx build --push --platform linux/amd64 --tag francesccoll/powerbacgwas:amd64 .
```
