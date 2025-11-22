CC = gcc
CFLAGS = -O2 -std=c11 -Wall -Wextra
TARGET = algdok

all: $(TARGET)

$(TARGET): algdok.c ../third_party/cjson/cJSON.c
	$(CC) $(CFLAGS) -I../third_party/cjson -o $(TARGET) algdok.c ../third_party/cjson/cJSON.c -lm

clean:
	rm -f $(TARGET)
