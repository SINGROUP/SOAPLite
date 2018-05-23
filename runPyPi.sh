sudo rm -r build/
sudo rm -r dist/
sudo python setup.py bdist_wheel
sudo python setup.py sdist
twine upload dist/*
