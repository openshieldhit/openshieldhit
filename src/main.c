#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[]) {
    /* Handle --version flag */
    if (argc > 1 && (strcmp(argv[1], "--version") == 0 || strcmp(argv[1], "-v") == 0)) {
        printf("OpenShieldHIT version %s\n", OSH_VERSION);
        return 0;
    }

    /* Handle --help flag */
    if (argc > 1 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0)) {
        printf("Usage: %s [OPTIONS]\n", argv[0]);
        printf("OpenShieldHIT - Monte Carlo Particle Transport\n\n");
        printf("OPTIONS:\n");
        printf("  --version, -v     Print version information\n");
        printf("  --help, -h        Show this help message\n");
        return 0;
    }

    printf("Hello, World!\n");
    return 0;
}