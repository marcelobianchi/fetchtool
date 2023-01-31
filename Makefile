all:
	echo "Choose one:"
	echo ""
	echo "make pipinstall"
	echo "make pipremove"
	echo ""

pipinstall:
	pip3 install --no-deps -e ./

pipremove:
	pip3 uninstall --yes fetchtool

