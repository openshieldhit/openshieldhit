#ifndef OSH_FILE_H
#define OSH_FILE_H

#include <stddef.h>
#include <stdio.h>

/**
 * @brief Internal structure to track an open text file with metadata
 */
struct oshfile {
    FILE *fp;
    long int *map;
    char *filename;
    size_t map_size;
    int map_len;
    int lineno;
};

/**
 * @brief Open a file for reading, allocate memory and initialize metadata.
 *
 * @details This function opens a file in read mode and initializes an oshfile_t
 * structure.
 *
 * @param filename Name of the file to open
 * @returns Pointer to an oshfile_t structure containing file metadata
 *
 * @author Niels Bassler
 */
struct oshfile *osh_fopen(char const *filename);

/**
 * @brief Close the file and clean up internal state

 * @details This function closes the file associated with the oshfile_t
 structure,
 *         frees the allocated memory for the filename, and deallocates the
 oshfile_t structure
 *         itself. It is safe to call this function with a NULL pointer.
 *
 * @param oshf Pointer to the oshfile_t structure to close
 *
 * @author Niels Bassler
 *
 */
void osh_fclose(struct oshfile *oshf);

/**
 * @brief Get the current line number in the file.
 *
 * @details This function returns the current line number based on the file
 * pointer position. It uses a binary search on the map of line offsets to find
 * the correct line number.
 *
 * @param oshf Pointer to the oshfile_t structure
 * @returns The current line number, or -1 if an error occurs
 *
 * @author Niels Bassler
 */
int osh_file_lineno(const struct oshfile *oshf);

/* Normalize path separators to '/' in-place.
 *
 * On Windows, the Win32 API accepts both '/' and '\' so normalizing to '/'
 * once in the thin client (main.c) lets all library code assume forward slashes
 * without per-subsystem #ifdef blocks.
 *
 * On non-Windows platforms this is a no-op. Safe to call with NULL. */
void osh_path_normalize(char *path);

/* Extract the directory portion of a path into a newly allocated string.
 *
 * Returns a malloc'd copy of everything up to (but not including) the last '/'.
 * Returns NULL if path is NULL, contains no separator, or allocation fails.
 * Caller owns the returned string and must free() it. */
char *osh_path_dirname(char const *path);

#endif /* OSH_FILE_H */