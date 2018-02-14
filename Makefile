.PHONY:	all clean allclean

all: | output
	@echo "\n*********** Build information ***********"
	@echo "  Debug is set to: [$(DEBUG)],"
	@echo "  Set it to 1 to enable a debug build."
	@echo "  For example: make clean; make DEBUG=1"
	@echo "*****************************************\n"
	$(MAKE) -C src

clean:
	$(MAKE) -C src clean

allclean:
	$(MAKE) -C src allclean

# Ensure the output folder exists
output:
		@mkdir -p output
