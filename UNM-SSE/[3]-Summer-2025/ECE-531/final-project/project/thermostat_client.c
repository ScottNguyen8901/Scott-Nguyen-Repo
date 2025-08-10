// thermostat_client.c
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#include <getopt.h>
#include <curl/curl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <fcntl.h>
#include <errno.h>

#define MAX_SETPOINTS 24
#define LOOP_DELAY 5

// Configuration structure
typedef struct {
    char temp_file[256];
    char status_file[256];
    char heater_log[256];
    char schedule_file[256];
    char server_url[256];
} Config;

typedef struct {
    int minutes;
    float temp;
} Setpoint;

static Setpoint schedule_[MAX_SETPOINTS];
static int schedule_count = 0;
static bool run_as_daemon = false;
static time_t sched_mtime = 0;

// ------------------------- Utilities -------------------------

static int clampi(int v, int lo, int hi) { return (v < lo) ? lo : (v > hi) ? hi : v; }

// Convert "HH:MM" to minutes since midnight
static int parse_time_to_minutes(const char* time_str) {
    int h = 0, m = 0;
    if (!time_str || sscanf(time_str, "%d:%d", &h, &m) != 2) return 0;
    h = clampi(h, 0, 23);
    m = clampi(m, 0, 59);
    return h * 60 + m;
}

// Insertion sort by minutes (small N)
static void sort_schedule(Setpoint* arr, int n) {
    for (int i = 1; i < n; i++) {
        Setpoint key = arr[i];
        int j = i - 1;
        while (j >= 0 && arr[j].minutes > key.minutes) {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = key;
    }
}

// ------------------------- Schedule -------------------------

// Load schedule file
static void load_schedule(const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Failed to open schedule '%s': %s\n", filename, strerror(errno));
        return;
    }

    int count = 0;
    char time_str[32];
    float temp;
    while (count < MAX_SETPOINTS && fscanf(fp, "%31s %f", time_str, &temp) == 2) {
        schedule_[count].minutes = parse_time_to_minutes(time_str);
        schedule_[count].temp = temp;
        count++;
    }
    fclose(fp);

    sort_schedule(schedule_, count);
    schedule_count = count;
}

// Hot-reload if schedule file mtime changed
static void maybe_reload_schedule(const char* path) {
    struct stat st;
    if (stat(path, &st) == 0) {
        if (st.st_mtime != sched_mtime) {
            load_schedule(path);
            sched_mtime = st.st_mtime;
        }
    } else {

    }
}

// Find current setpoint based on system time
static float get_current_setpoint(void) {
    // Default if no schedule loaded yet
    if (schedule_count <= 0) return 70.0f;

    time_t now = time(NULL);
    struct tm* t = localtime(&now);
    int current_minutes = t->tm_hour * 60 + t->tm_min;

    float latest_temp = schedule_[0].temp;
    for (int i = 0; i < schedule_count; i++) {
        if (current_minutes >= schedule_[i].minutes) {
            latest_temp = schedule_[i].temp;
        } else {
            break;
        }
    }
    return latest_temp;
}

// ------------------------- Config -------------------------

static void load_config(const char* config_path, Config* cfg) {
    FILE* fp = fopen(config_path, "r");
    if (!fp) {
        fprintf(stderr, "Failed to open config file '%s': %s\n", config_path, strerror(errno));
        exit(1);
    }

    while (!feof(fp)) {
        char key[256], value[256];
        if (fscanf(fp, "%255s %255s", key, value) != 2) continue;

        if (strcmp(key, "temperature_file") == 0)
            strncpy(cfg->temp_file, value, sizeof(cfg->temp_file));
        else if (strcmp(key, "status_file") == 0)
            strncpy(cfg->status_file, value, sizeof(cfg->status_file));
        else if (strcmp(key, "heater_log") == 0)
            strncpy(cfg->heater_log, value, sizeof(cfg->heater_log));
        else if (strcmp(key, "schedule_file") == 0)
            strncpy(cfg->schedule_file, value, sizeof(cfg->schedule_file));
        else if (strcmp(key, "server_endpoint") == 0)
            strncpy(cfg->server_url, value, sizeof(cfg->server_url));
    }

    fclose(fp);

    cfg->temp_file[sizeof(cfg->temp_file) - 1] = '\0';
    cfg->status_file[sizeof(cfg->status_file) - 1] = '\0';
    cfg->heater_log[sizeof(cfg->heater_log) - 1] = '\0';
    cfg->schedule_file[sizeof(cfg->schedule_file) - 1] = '\0';
    cfg->server_url[sizeof(cfg->server_url) - 1] = '\0';
}

// ------------------------- HTTP -------------------------

static void report_to_server(const char* url, float temp, const char* heater_state) {
    if (!url || url[0] == '\0') return;

    CURL* curl = curl_easy_init();
    if (!curl) return;

    char json[128];
    snprintf(json, sizeof(json), "{\"temperature\": %.2f, \"heater\": \"%s\"}", temp, heater_state);

    struct curl_slist* headers = NULL;
    headers = curl_slist_append(headers, "Content-Type: application/json");

    curl_easy_setopt(curl, CURLOPT_URL, url);
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, json);
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);

    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 3L);
    curl_easy_setopt(curl, CURLOPT_CONNECTTIMEOUT, 2L);
    curl_easy_setopt(curl, CURLOPT_NOSIGNAL, 1L);

    CURLcode res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
        fprintf(stderr, "HTTP POST failed: %s\n", curl_easy_strerror(res));
    }

    curl_slist_free_all(headers);
    curl_easy_cleanup(curl);
}

// ------------------------- Daemon -------------------------

static void daemonize(void) {
    pid_t pid = fork();
    if (pid < 0) exit(EXIT_FAILURE);
    if (pid > 0) exit(EXIT_SUCCESS);

    if (setsid() < 0) exit(EXIT_FAILURE);

    pid = fork();
    if (pid < 0) exit(EXIT_FAILURE);
    if (pid > 0) exit(EXIT_SUCCESS);

    umask(0);

    if (chdir("/") != 0) exit(EXIT_FAILURE);

    struct rlimit rl;
    if (getrlimit(RLIMIT_NOFILE, &rl) != 0) {
        rl.rlim_max = 1024;
    }
    for (int fd = 0; fd < (int)rl.rlim_max; fd++) {
        close(fd);
    }

    int fdnull = open("/dev/null", O_RDWR);
    if (fdnull < 0) exit(EXIT_FAILURE);
    if (dup2(fdnull, STDIN_FILENO)  < 0) exit(EXIT_FAILURE);
    if (dup2(fdnull, STDOUT_FILENO) < 0) exit(EXIT_FAILURE);
    if (dup2(fdnull, STDERR_FILENO) < 0) exit(EXIT_FAILURE);
    if (fdnull > STDERR_FILENO) close(fdnull);
}

// ------------------------- UI -------------------------

static void print_help(void) {
    printf("Usage: ./thermostat_client [-c config_file] [-d] [--help]\n");
    printf("  -c, --config_file <file> : Path to configuration file (default: config.cfg)\n");
    printf("  -d, --daemon             : Run as background daemon\n");
    printf("  -h, --help               : Show this help message\n");
}

// ------------------------- Main -------------------------

int main(int argc, char* argv[]) {
    Config cfg = {0};
    // Defaults — override in config.cfg
    strncpy(cfg.temp_file, "/var/log/temperature", sizeof(cfg.temp_file));
    strncpy(cfg.status_file, "/tmp/status", sizeof(cfg.status_file));
    strncpy(cfg.heater_log, "/var/log/heater", sizeof(cfg.heater_log));
    strncpy(cfg.schedule_file, "schedule.txt", sizeof(cfg.schedule_file));
    strncpy(cfg.server_url, "http://localhost:8000/status", sizeof(cfg.server_url));

    const char* config_file = "config.cfg";

    static struct option long_options[] = {
        {"config_file", required_argument, 0, 'c'},
        {"daemon",      no_argument,       0, 'd'},
        {"help",        no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "c:dh", long_options, NULL)) != -1) {
        switch (opt) {
            case 'c': config_file = optarg; break;
            case 'd': run_as_daemon = true; break;
            case 'h': print_help(); return 0;
            default:  print_help(); return 1;
        }
    }

    if (run_as_daemon) {
        daemonize();
    }

    load_config(config_file, &cfg);

    maybe_reload_schedule(cfg.schedule_file);

    float current_temp = 0.0f;
    char last_action[4] = "xxx";

    while (1) {
        maybe_reload_schedule(cfg.schedule_file);

        FILE* temp_fp = fopen(cfg.temp_file, "r");
        if (!temp_fp || fscanf(temp_fp, "%f", &current_temp) != 1) {
            fprintf(stderr, "Failed to read temperature from '%s': %s\n", cfg.temp_file, strerror(errno));
            if (temp_fp) fclose(temp_fp);
            sleep(LOOP_DELAY);
            continue;
        }
        fclose(temp_fp);

        float setpoint = get_current_setpoint();
        const char* action = (current_temp < setpoint) ? "on" : "off";

        if (strcmp(last_action, action) != 0) {
            time_t now = time(NULL);

            FILE* log_fp = fopen(cfg.heater_log, "a");
            if (log_fp) {
                fprintf(log_fp, "%s : %ld\n", action, (long)now);
                fclose(log_fp);
            } else {
                fprintf(stderr, "Failed to write heater_log '%s': %s\n", cfg.heater_log, strerror(errno));
            }

            FILE* status_fp = fopen(cfg.status_file, "w");
            if (status_fp) {
                fprintf(status_fp, "%s\n", action);
                fclose(status_fp);
            } else {
                fprintf(stderr, "Failed to write status_file '%s': %s\n", cfg.status_file, strerror(errno));
            }

            report_to_server(cfg.server_url, current_temp, action);
            strncpy(last_action, action, sizeof(last_action));

            if (!run_as_daemon) {
                printf("[%.1f °C] Setpoint %.1f → Heater turned %s at %ld\n",
                       current_temp, setpoint, action, (long)now);
            }
        }

        report_to_server(cfg.server_url, current_temp, action);
        sleep(LOOP_DELAY);
    }

    return 0;
}
