CC = gcc
CFLAGS = -O2 -std=c11 -Wall -Wextra
TARGET = graph_to_json

all: $(TARGET)

$(TARGET): graph_to_json.c ../third_party/cjson/cJSON.c
	$(CC) $(CFLAGS) -I../third_party/cjson -o $(TARGET) graph_to_json.c ../third_party/cjson/cJSON.c

clean:
	rm -f $(TARGET)

# Note: system package names: 'libcjson-dev' (Debian/Ubuntu). Alternatively compile cJSON.c into project.
