TARGET=fagin
SOURCE=src
CFLAGS=-Wall -g -o $(TARGET)
CC=cc

# If anything changes, recompile the whole thing
all:
	$(CC) $(CFLAGS) $(SOURCE)/*c $(SOURCE)/*h || rm -f $(SOURCE)/*.gch

clean:
	rm -f $(TARGET) $(SOURCE)/*.gch
