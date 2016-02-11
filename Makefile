TARGET=main
SOURCE=src
CFLAGS=-Wall -g -o $(TARGET)
CC=cc

# If anything changes, recompile the whole thing
all:
	$(CC) $(CFLAGS) $(SOURCE)/*.c $(SOURCE)/itree/*.c $(SOURCE)/*.h $(SOURCE)/itree/*.h || rm -f $(SOURCE)/**.gch

docs:
	doxygen Doxyfile

clean:
	rm -f $(TARGET) $(SOURCE)/*.gch
