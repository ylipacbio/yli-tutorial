all: meson nj execute

meson:
	meson  --default-library=static build-static .
nj:
	ninja -C build-static -v -j 20 |tee log
	tail log
execute:
	./build-static/my
