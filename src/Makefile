CC = gcc

DEPS_PREFIX ?= /bi/software/static
HTSLIB_DIR ?= $(DEPS_PREFIX)
HTSLIB_INCLUDE ?= $(HTSLIB_DIR)
HTSLIB_LIB ?= $(HTSLIB_DIR)
BUILD_DIR ?= build
TARGET ?= $(BUILD_DIR)/cluster_identifier
SRC = cluster_identifier.c

# Set STATIC=0 for a normal dynamic build.
STATIC ?= 1

# Match these to the htslib build:
#   USE_CURL=1       if htslib was built with libcurl support
#   USE_LIBDEFLATE=1 if htslib was built with libdeflate support
USE_CURL ?= 0
USE_LIBDEFLATE ?= 0

CPPFLAGS += -I$(HTSLIB_INCLUDE) -I$(DEPS_PREFIX)/include
CFLAGS += -Wall -O2
LDFLAGS += -L$(HTSLIB_LIB) -L$(DEPS_PREFIX)/lib

ifeq ($(STATIC),1)
  LDFLAGS += -static
endif

LDLIBS = -lhts

ifeq ($(USE_LIBDEFLATE),1)
  LDLIBS += -ldeflate
endif

LDLIBS += -lz -llzma -lbz2

ifeq ($(USE_CURL),1)
  LDLIBS += -lcurl
endif

LDLIBS += -lcrypto -lpthread -lm -ldl

.PHONY: all clean

all: $(TARGET)
	@echo $(TARGET)

$(TARGET): $(SRC)
	@mkdir -p $(@D)
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) $(TARGET)
