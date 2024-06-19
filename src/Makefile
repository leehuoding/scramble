## Compilation options/libs
CC = gcc
CFLAGS = -Wall -O2 -I/path/to/libdeflate -I/path/to/htslib
LIBS = -L/path/to/libdeflate/build -L/path/to/htslib -lpthread -llzma -lbz2 -lcurl -lcrypto -lm -lz -lhts -ldeflate 

## source and exe names
SRC = cluster_identifier.c
BASEMAIN = cluster_identifier
BUILD_DIR = $(dir $(realpath $(firstword $(MAKEFILE_LIST))))
MAIN = $(BUILD_DIR)/$(BASEMAIN)

.PHONY: clean

all: $(MAIN)
	@echo $(MAIN)

$(MAIN):
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(SRC) $(LIBS)

clean:
	$(RM) $(BUILD_DIR)/$(BASEMAIN)
