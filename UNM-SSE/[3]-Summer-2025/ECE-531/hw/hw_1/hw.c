// cc -o hw hw.c
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <curl/curl.h>

#define OK 0
#define INIT_ERR 1
#define REQ_ERR 2
#define ARG_ERR 3

void print_help() {
    printf("Usage:\n");
    printf("  hw [OPTIONS] MESSAGE\n\n");
    printf("Options:\n");
    printf("  -u, --url URL           Set the target URL (e.g., http://localhost:8000)\n");
    printf("  -g, --get               Perform HTTP GET request\n");
    printf("  -o, --post              Perform HTTP POST request\n");
    printf("  -p, --put               Perform HTTP PUT request\n");
    printf("  -d, --delete            Perform HTTP DELETE request\n");
    printf("  -h, --help              Display this help message\n");
}

size_t write_callback(void *contents, size_t size, size_t nmemb, void *userp) {
    fwrite(contents, size, nmemb, stdout);
    return size * nmemb;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        print_help();
        return ARG_ERR;
    }

    char *url = NULL;
    char *method = NULL;
    int i = 1;

    // Gather method and URL
    while (i < argc) {
        if ((strcmp(argv[i], "--url") == 0 || strcmp(argv[i], "-u") == 0) && i + 1 < argc) {
            url = argv[++i];
        } else if (strcmp(argv[i], "--get") == 0 || strcmp(argv[i], "-g") == 0) {
            method = "GET";
        } else if (strcmp(argv[i], "--post") == 0 || strcmp(argv[i], "-o") == 0) {
            method = "POST";
        } else if (strcmp(argv[i], "--put") == 0 || strcmp(argv[i], "-p") == 0) {
            method = "PUT";
        } else if (strcmp(argv[i], "--delete") == 0 || strcmp(argv[i], "-d") == 0) {
            method = "DELETE";
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            print_help();
            return OK;
        } else {
            break;
        }
        i++;
    }

    // Remaining args = message (if needed)
    char data[2048] = {0};
    while (i < argc) {
        strcat(data, argv[i++]);
        if (i < argc) strcat(data, " ");
    }

    if (!url || !method) {
        fprintf(stderr, "Error: URL and method required.\n");
        print_help();
        return ARG_ERR;
    }

    CURL *curl = curl_easy_init();
    if (!curl) return INIT_ERR;

    CURLcode res;
    curl_easy_setopt(curl, CURLOPT_URL, url);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);

    if (strcmp(method, "POST") == 0) {
        if (strlen(data) == 0) return ARG_ERR;
        curl_easy_setopt(curl, CURLOPT_POST, 1L);
        curl_easy_setopt(curl, CURLOPT_POSTFIELDS, data);
    } else if (strcmp(method, "PUT") == 0) {
        if (strlen(data) == 0) return ARG_ERR;
        curl_easy_setopt(curl, CURLOPT_CUSTOMREQUEST, "PUT");
        curl_easy_setopt(curl, CURLOPT_POSTFIELDS, data);
    } else if (strcmp(method, "DELETE") == 0) {
        curl_easy_setopt(curl, CURLOPT_CUSTOMREQUEST, "DELETE");
        if (strlen(data) > 0) {
        // Removed request to body in DELETE request. Added print statement for previous misunderstanding on my part
        fprintf(stderr, "Warning: DELETE requests typically do not send a body. Ignoring message.\n");
        }
    }

    res = curl_easy_perform(curl);

    long response_code = 0;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &response_code);

    printf("\nHTTP Status: %ld\n", response_code);

    curl_easy_cleanup(curl);
    return (res == CURLE_OK) ? OK : REQ_ERR;
}