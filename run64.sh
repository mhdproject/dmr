tar czvf ath.tgz athena-graphics-3d
tar czvf par.tgz par
docker build --file Dockerfile.64 --tag mhd64 .
docker run -it mhd64 /bin/bash
