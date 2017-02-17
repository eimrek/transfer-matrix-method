
ifeq ($(OS), Windows_NT)
    	# Windows environment
	APP = tmm.exe
	XCFLAGS = 
	XLFLAGS = -lmingw32
else
	# Linux environment
	APP = tmm
	XCFLAGS =
	XLFLAGS = -lm
endif

CC = g++

CFLAGS = -g -Wall --std=c++11 $(XCFLAGS)
LFLAGS = $(XLFLAGS)

$(APP): main.cpp
	@mkdir -p data
	$(CC) $(CFLAGS) $(LFLAGS) $^ -o $@

clean:
	-rm -f *.o
	-rm -f $(APP)
