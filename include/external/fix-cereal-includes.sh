cd cereal
for f in *; do [[ -d "$f" ]] || sed -i '' -e 's/#include <cereal\/\([^>]*\)>/#include "\1"/g' "$f"; done
for subdirs in $(\ls -d */); do
    cd ${subdirs%%/}
    for f in *; do [[ -d "$f" ]] || sed -i '' -e 's/#include <cereal\/\([^>]*\)>/#include "..\/\1"/g' "$f"; done
    for subdirs in $(\ls -d */); do
	cd ${subdirs%%/}
	for f in *; do [[ -d "$f" ]] || sed -i '' -e 's/#include <cereal\/\([^>]*\)>/#include "..\/..\/\1"/g' "$f"; done
	cd ..
    done
    cd ..
done
cd ..
