/*Tugas Besar EL2208 Praktikum Pemecahan Masalah dengan C 2023/2024
*Kelompok         : F3
*Anggota          : Muhammad Zaki Fazansyah  (18322018) 
                    Andika Rama Ardhi Nugraha  (18322019) 
                    Yazid Nur Hidayat         (18322020) 
                    Requint Nurliawan         (18322021) 
                    Muhammad Rafi Ar-Rantisi (18322022) 
                    Mayesa Adhisty          (18322023) 
                    Musthafa Ibrahim       (18322024) 
*Asisten (NIM)    : Muhammad Daris Nurhakim (13220047)
*Nama File        : Tubes PPMC_F3.c
*Deskripsi        : The Travelling Salesman Problem : menentukan rute terbaik dan jarak terdekat dengan melewati semua kota
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <sys/time.h>

#define EARTH_RADIUS 6371.0
#define PI 3.14159265358979323846

// Greedy define
#define MAKS_KOTA 15

// ACO define
#define MAX_CITY 1000
#define MAX_NAME 100
#define MAX_ANT 100
#define MAX_ITERATIONS 100
#define ALPHA 3.0
#define BETA 4.0
#define Q 10000.0
#define RHO 0.00001

// Genetics define
#define POP_SIZE 100
#define NUM_CITIES_init 8
#define MUTATION_RATE 0.1
#define MAX_GEN 1000


// START OF DFS
typedef struct {
    char name[50];
    double lat;
    double lon;
} Kota;

Kota daftar_kota[15];
int n;
double jarak_min = DBL_MAX;
int best_route[15];
int current_route[15];
int visited[15];

double to_radians_dfs(double degree) {
    return degree * (PI / 180);
}

double haversine_dfs(double lat1, double lon1, double lat2, double lon2) {
    double dlat = to_radians_dfs(lat2 - lat1);
    double dlon = to_radians_dfs(lon2 - lon1);
    double a = sin(dlat/2) * sin(dlat/2) + cos(to_radians_dfs(lat1)) * cos(to_radians_dfs(lat2)) * sin(dlon/2) * sin(dlon/2);
    double c = 2 * asin(sqrt(a));
    return EARTH_RADIUS * c;
}

void jarak_terpendek_dfs(int current, int count, double jarak) {
    if (count == n) {
        jarak += haversine_dfs(daftar_kota[current].lat, daftar_kota[current].lon, daftar_kota[current_route[0]].lat, daftar_kota[current_route[0]].lon);
        if (jarak < jarak_min) {
            jarak_min = jarak;
            memcpy(best_route, current_route, sizeof(current_route));
        }
        return;
    }
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            visited[i] = 1;
            current_route[count] = i;
            jarak_terpendek_dfs(i, count + 1, jarak + haversine_dfs(daftar_kota[current].lat, daftar_kota[current].lon, daftar_kota[i].lat, daftar_kota[i].lon));
            visited[i] = 0;
        }
    }
}

int find_city_index_dfs(char* start_city) {
    for (int i = 0; i < n; i++) {
        if (strcmp(daftar_kota[i].name, start_city) == 0) {
            return i;
        }
    }
    return -1;
}

int dfs(){
    FILE *file;
    char filename[100], start_city[50];
    
    printf("Enter list of cities file name: ");
    scanf("%s", filename);
    printf("Enter starting point: ");
    scanf("%s", start_city);
    clock_t start_time, end_time;

    file = fopen(filename, "r");
    if (!file) {
        printf("Unable to open file.");
        return 1;
    }

    n = 0;
    while (fscanf(file, "%[^,],%lf,%lf\n", daftar_kota[n].name, &daftar_kota[n].lat, &daftar_kota[n].lon) == 3) {
        n++;
    }
    fclose(file);

    int start_index = find_city_index_dfs(start_city);
    if (start_index == -1) {
        printf("Starting city not found in the list.\n");
        return 1;
    }

    memset(visited, 0, sizeof(visited));
    current_route[0] = start_index;
    visited[start_index] = 1;

    start_time = clock();
    jarak_terpendek_dfs(start_index, 1, 0);
    end_time = clock();

    printf("The shortest route is:\n");
    for (int i = 0; i < n; i++) {
        printf("%s -> ", daftar_kota[best_route[i]].name);
    }
    printf("%s\n", daftar_kota[start_index].name);
    printf("Best Route_BFS Distance: %f km\n", jarak_min);
    printf("Time elapsed: %f s\n", ((double)(end_time - start_time)) / CLOCKS_PER_SEC);

    return 0;
}
// END OF DFS

// START OF BFS
typedef struct {
    char name[50];  //Nama Kota 
    double latitude; //Garis lintang kota (double)
    double longitude; //Garis bujur kota (double)
} City_BFS; //Menyimpan informasi ttg kota

typedef struct {
    int *path; //Menyimpan urutan indeks kota yang dikunjungi
    int length; //Panjang dr path (jumlah kota yang dikunjungi)
    double distance; //Jarak total dari rute yang dilewati
} Route_BFS; //Menyimpan informasi ttg rute 

double toRadians_bfs(double degree) { //mengkonversi derajat ke radian 
    return degree * (PI / 180.0);
}

double haversineDistance_bfs(double lat1, double lon1, double lat2, double lon2) { //Menghitung jarak dua titik dalam derajat (Rumus Haversine)
    double dLat = toRadians_bfs(lat2 - lat1);
    double dLon = toRadians_bfs(lon1 - lon2);
    double a = sin(dLat / 2) * sin(dLat / 2) +
               cos(toRadians_bfs(lat1)) * cos(toRadians_bfs(lat2)) *
               sin(dLon / 2) * sin(dLon / 2);
    double c = 2 * asin(sqrt(a));
    return EARTH_RADIUS * c;
}

double calculateTotalDistance_bfs(City_BFS *cities, int *path, int n) { //Menghitung jarak total dari rute yang diberikan
    double totalDistance = 0.0;
    for (int i = 0; i < n - 1; ++i) {
        totalDistance += haversineDistance_bfs(cities[path[i]].latitude, cities[path[i]].longitude, 
                                           cities[path[i + 1]].latitude, cities[path[i + 1]].longitude);
    }
    // Add distance back to the starting city
    totalDistance += haversineDistance_bfs(cities[path[n - 1]].latitude, cities[path[n - 1]].longitude, 
                                       cities[path[0]].latitude, cities[path[0]].longitude);
    return totalDistance;
}

void bfsTSP(City_BFS *cities, int n, int startIdx, Route_BFS *bestRoute) { //Algoritma BFS 
    int queueCapacity = 1000000;
    int **queue = malloc(queueCapacity * sizeof(int*));
    int front = 0, rear = 0;
    
    queue[rear] = malloc((n + 1) * sizeof(int));
    queue[rear][0] = startIdx;
    queue[rear][n] = 1;  // path length
    rear++;
    
    bestRoute->distance = INFINITY;

    while (front < rear) {
        int *currentPath = queue[front];
        int pathLength = currentPath[n];
        front++;
        
        if (pathLength == n) {
            double currentDistance = calculateTotalDistance_bfs(cities, currentPath, n);
            if (currentDistance < bestRoute->distance) {
                bestRoute->distance = currentDistance;
                memcpy(bestRoute->path, currentPath, n * sizeof(int));
                bestRoute->length = n;
            }
            free(currentPath);
            continue;
        }

        for (int i = 0; i < n; ++i) {
            int alreadyVisited = 0;
            for (int j = 0; j < pathLength; ++j) {
                if (currentPath[j] == i) {
                    alreadyVisited = 1;
                    break;
                }
            }
            if (!alreadyVisited) {
                queue[rear] = malloc((n + 1) * sizeof(int));
                memcpy(queue[rear], currentPath, pathLength * sizeof(int));
                queue[rear][pathLength] = i;
                queue[rear][n] = pathLength + 1;
                rear++;
            }
        }
        free(currentPath);
    }
    for (int i = front; i < rear; ++i) {
        free(queue[i]);
    }
    free(queue);
}

int bfs() {
    char filename[100], startCity[50]; 
    //filename : menyimpan nama file yang berisi daftar kota
    //startCity : menyimpan nama kota awal 
    printf("Enter list of cities file name: ");
    scanf("%s", filename);
    getchar(); // Consume newline left by scanf
    printf("Enter starting point: ");
    fgets(startCity, sizeof(startCity), stdin);
    startCity[strcspn(startCity, "\n")] = 0; // Remove newline character

    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error opening file!\n");
        return 1;
    }

    City_BFS *cities = malloc(100 * sizeof(City_BFS));
    int cityCount = 0; //cityCount : menyimpan jumlah kota yang telah dibaca dari file 
    int cityCapacity = 100; //cityCapacity : menyimpan kapasitas saat ini dari array cities 

    while (fscanf(file, "%49[^,],%lf,%lf\n", cities[cityCount].name, &cities[cityCount].latitude, &cities[cityCount].longitude) != EOF) {
        cityCount++;
        if (cityCount >= cityCapacity) {
            cityCapacity *= 2;
            cities = realloc(cities, cityCapacity * sizeof(City_BFS));
        }
    }
    fclose(file);

    int startIdx = -1;
    for (int i = 0; i < cityCount; ++i) {
        if (strcmp(cities[i].name, startCity) == 0) {
            startIdx = i;
            break;
        }
    }

    if (startIdx == -1) {
        printf("Starting city not found!\n");
        free(cities);
        return 1;
    }

    Route_BFS bestRoute; //menyimpan rute terbaik yang ditemukan
    bestRoute.path = malloc(cityCount * sizeof(int));
    bestRoute.length = 0;
    bestRoute.distance = INFINITY;

    clock_t begin = clock(); //mengukur waktu yang dibutuhkan
    bfsTSP(cities, cityCount, startIdx, &bestRoute);
    clock_t end = clock(); //mengukur waktu yang dibutuhkan
    double timeSpent = (double)(end - begin) / CLOCKS_PER_SEC; //menyimpan waktu yang dibutuhkan

    printf("Best route found:\n");
    printf("%s", cities[startIdx].name);
    for (int i = 0; i < bestRoute.length; ++i) {
        if (bestRoute.path[i] != startIdx) {
            printf(" -> %s", cities[bestRoute.path[i]].name);
        }
    }
    printf(" -> %s\n", cities[startIdx].name); // Print start city at the end

    printf("Best route distance: %lf km\n", bestRoute.distance);
    printf("Time elapsed: %lf s\n", timeSpent);

    free(bestRoute.path);
    free(cities);

    return 0;
}
//END OF BFS

//START OF GREEDY

typedef struct {
    char nama[50];
    double latitude;
    double longitude;
} Kotahehe; // Ubah nama struct menjadi Kotahehe

double haversine(double lat1, double lon1, double lat2, double lon2) {
    double phi1 = lat1 * PI / 180.0;
    double phi2 = lat2 * PI / 180.0;
    double delta_phi = (lat2 - lat1) * PI / 180.0;
    double delta_lambda = (lon2 - lon1) * PI / 180.0;

    double a = sin(delta_phi / 2.0) * sin(delta_phi / 2.0) +
               cos(phi1) * cos(phi2) * sin(delta_lambda / 2.0) * sin(delta_lambda / 2.0);
    double c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a));

    return EARTH_RADIUS * c;
}

int cari_kota_terdekat(int kota_saat_ini, int n, int dikunjungi[], Kotahehe kota[]) {
    double jarak_min = DBL_MAX;
    int kota_terdekat = -1;

    for (int i = 0; i < n; i++) {
        if (!dikunjungi[i]) {
            double jarak = haversine(kota[kota_saat_ini].latitude, kota[kota_saat_ini].longitude,
                                     kota[i].latitude, kota[i].longitude);
            if (jarak < jarak_min) {
                jarak_min = jarak;
                kota_terdekat = i;
            }
        }
    }

    return kota_terdekat;
}

void tsp_greedy(int kota_awal, int n, Kotahehe kota[]) {
    int dikunjungi[n];
    for (int i = 0; i < n; i++) {
        dikunjungi[i] = 0;
    }

    int kota_saat_ini = kota_awal;
    dikunjungi[kota_saat_ini] = 1;

    double jarak_total = 0;
    printf("Best route found:\n %s", kota[kota_saat_ini].nama);

    for (int i = 1; i < n; i++) {
        int kota_berikutnya = cari_kota_terdekat(kota_saat_ini, n, dikunjungi, kota);
        jarak_total += haversine(kota[kota_saat_ini].latitude, kota[kota_saat_ini].longitude,
                                 kota[kota_berikutnya].latitude, kota[kota_berikutnya].longitude);
        kota_saat_ini = kota_berikutnya;
        dikunjungi[kota_saat_ini] = 1;
        printf(" -> %s", kota[kota_saat_ini].nama);
    }

    jarak_total += haversine(kota[kota_saat_ini].latitude, kota[kota_saat_ini].longitude,
                             kota[kota_awal].latitude, kota[kota_awal].longitude);
    printf(" -> %s\n", kota[kota_awal].nama);

    printf("Jarak total: %.5f km\n", jarak_total);
}

int baca_kota_dari_csv(const char *nama_file, Kotahehe kota[]) {
    FILE *file = fopen(nama_file, "r");
    if (!file) {
        printf("Error opening file! %s\n", nama_file);
        return -1;
    }

    char baris[100];
    int count = 0;

    while (fgets(baris, sizeof(baris), file)) {
        if (sscanf(baris, "%49[^,],%lf,%lf", kota[count].nama, &kota[count].latitude, &kota[count].longitude) == 3) {
            count++;
            if (count >= MAKS_KOTA) {
                printf("List exceed the maximum number (%d) of cities.\n", MAKS_KOTA);
                break;
            }
        }
    }

    fclose(file);
    return count;
}

int greedy() {
    char nama_file[100];
    printf("Masukkan nama file CSV: ");
    scanf("%99s", nama_file);

    Kotahehe kota[MAKS_KOTA]; // Ubah tipe data struct menjadi Kotahehe
    int n = baca_kota_dari_csv(nama_file, kota);
    if (n <= 0) {
        return 1;
    }

    char nama_kota_awal[50];
    printf("Enter starting point: ");
    getchar();
    fgets(nama_kota_awal, sizeof(nama_kota_awal), stdin);
    nama_kota_awal[strcspn(nama_kota_awal, "\n")] = '\0';

    int kota_awal = -1;
    for (int i = 0; i < n; i++) {
        if (strcmp(kota[i].nama, nama_kota_awal) == 0) {
            kota_awal = i;
            break;
        }
    }

    if (kota_awal == -1) {
        printf("City is not in the list.\n");
        return 1;
    }

    clock_t mulai = clock();

    tsp_greedy(kota_awal, n, kota);

    clock_t selesai = clock();
    double waktu_eksekusi = (double)(selesai - mulai) / CLOCKS_PER_SEC;

    printf("Time elapsed: %.6f s\n", waktu_eksekusi);

    return 0;
}


//END OF GREEDY

//START OF BRUTE FORCE
int size;
int vis[15], permutation[15], n;
int best_route[15];
float min_cost = 1e9;
char namaKota[15][255];
float matriks[15][15];

int open_init(char *namaFile, float matriks[15][15], char namaKota[15][255], int *n){
    *n = 0;
    // membuka dan cek file
    FILE *file = fopen(namaFile, "r");
    if (file == NULL){
        printf("File tidak ditemukan.");
        return 0;
    }
    
    char line[255];
    float latitude[15];
    float longitude[15];
    
    //membaca file
    while (fgets(line, 255, file)){
        //pengambilan data dari baris file
        strcpy(namaKota[*n], strtok(line, ","));
        latitude[*n] = atof(strtok(NULL, ","));
        longitude[*n] = atof(strtok(NULL, "\n"));
        
        *n += 1;
    }
    
    // cek kesesuaian jumlah kota (4 <= jumlah kota <= 15)
    if (*n < 3 || *n > 15){
        printf("\nJumlah kota tidak sesuai.");
        return 0;
    }
    // mengisi matriks yang berisi jarak dari suatu titik ke titik lainnya
    for (int i = 0; i < *n; i++){
        for (int j = 0; j < *n; j++){
            matriks[i][j] = 2 * EARTH_RADIUS * asin(sqrt(pow(sin((latitude[j]-latitude[i])*PI/360), 2) + cos(latitude[j]*PI/180) * cos(latitude[i]*PI/180) * pow(sin((longitude[j]-longitude[i])*PI/360), 2)));
        }
    }
    return 1;
}

void generate(int x, float cost) {
    if(x == n) { 
        cost += matriks[permutation[n - 1]][permutation[0]]; // Tambahkan jarak kembali ke kota awal
        if(cost < min_cost) { 
            min_cost = cost;
            for(int i = 0; i < n; ++i) {
                best_route[i] = permutation[i];
            }
        }
    } else {
        for(int i = 0; i < n; ++i) {
            if(vis[i]) continue;
            vis[i] = 1;
            permutation[x] = i;
            generate(x + 1, cost + matriks[permutation[x-1]][i]);
            vis[i] = 0;
        }
    }
}

void find_starting_point(char *point, int *index) {
    for(int i = 0; i < n; ++i) {
        if(strcmp(point, namaKota[i]) == 0) {
            *index = i;
            return;
        }
    }
    *index = -1;
}

void output_best_route(){
    puts("Best route found:");
    for(int i = 0; i < n; ++i) {
        printf("%s -> ", namaKota[best_route[i]]);
    }

    printf("%s\n", namaKota[best_route[0]]);
    printf("Best route distance: %.5f km\n", min_cost);
}

int brute_force() {
    char namaFile[255];
    n = 0;
 
    printf("Enter list of cities file name: ");
    scanf("%s", namaFile);

    if (open_init(namaFile, matriks, namaKota, &n) == 0){
        return 0;
    }

    size = n; 

    char startingCity[100];
    printf("Enter starting point: ");
    getchar(); 
    fgets(startingCity, sizeof(startingCity), stdin);
    startingCity[strcspn(startingCity, "\n")] = '\0'; 

    int startingIndex;
    find_starting_point(startingCity, &startingIndex);

    if (startingIndex == -1) {
        printf("Starting city not found\n");
        return 0;
    }

    vis[startingIndex] = 1; 
    permutation[0] = startingIndex;
    clock_t start_time = clock();
    generate(1,0);
    clock_t end_time = clock();
    output_best_route();

    double time_spent = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Time elapsed: %.10f s\n", time_spent);

    return 0;
}
//END OF BRUTE FORCEC

//START OF BRANCH N BOUND
typedef struct {
    char name[50];
    double latitude;
    double longitude;
} City_bnb;

int numCities = 0;
City_bnb cities[15];
double dist[15][15];
int bestRoute[15];
double bestDistance = DBL_MAX;

// Fungsi untuk mengubah derajat ke radian
double toRadians_bnb(double degree) {
    return degree * (PI / 180);
}

// Rumus Haversine
double haversine_bnb(double lat1, double lon1, double lat2, double lon2) {
    double dLat = toRadians_bnb(lat2 - lat1);
    double dLon = toRadians_bnb(lon2 - lon1);
    lat1 = toRadians_bnb(lat1);
    lat2 = toRadians_bnb(lat2);
    double a = sin(dLat / 2) * sin(dLat / 2) +
               sin(dLon / 2) * sin(dLon / 2) * cos(lat1) * cos(lat2);
    double c = 2 * asin(sqrt(a));
    return EARTH_RADIUS * c;
}

// Fungsi membaca nama kota di file
void citiesFile_bnb(char *fileName) {
    FILE *file = fopen(fileName, "r");
    if (file == NULL) {
        printf("File tidak ditemukan: %s\n", fileName);
        exit(1);
    }
    while (fscanf(file, "%[^,],%lf,%lf\n", cities[numCities].name, &cities[numCities].latitude, &cities[numCities].longitude) != EOF) {
        // printf("Read city: %s, %lf, %lf\n", cities[numCities].name, cities[numCities].latitude, cities[numCities].longitude); // Debugging
        numCities++;
    }
    fclose(file);
}

// Fungsi untuk menghitung jarak
void computeDistances() {
    for (int i = 0; i < numCities; i++) {
        for (int j = 0; j < numCities; j++) {
            if (i != j) {
                dist[i][j] = haversine_bnb(cities[i].latitude, cities[i].longitude, cities[j].latitude, cities[j].longitude);
            } else {
                dist[i][j] = 0;
            }
        }
    }
}

// Fungsi rekursif branch and bound
void branchbound(int currentCity, int visitedCities, double currentDistance, int route[], int depth) {
    if (depth == numCities) {
        double totalDistance = currentDistance + dist[currentCity][route[0]];
        if (totalDistance < bestDistance) {
            bestDistance = totalDistance;
            for (int i = 0; i < numCities; i++) {
                bestRoute[i] = route[i];
            }
        }
        return;
    }

    for (int i = 0; i < numCities; i++) {
        if (!(visitedCities & (1 << i))) {
            route[depth] = i;
            branchbound(i, visitedCities | (1 << i), currentDistance + dist[currentCity][i], route, depth + 1);
        }
    }
}

// Fungsi untuk mengetahui indeks kota pertama
int findStartingCity(char *startCity) {
    for (int i = 0; i < numCities; i++) {
        if (strcmp(cities[i].name, startCity) == 0) {
            return i;
        }
    }
    return -1;
}

// Fungsi untuk mencetak rute terpendek
void printRoute() {
    printf("Best route found:\n");
    for (int i = 0; i < numCities; i++) {
        printf("%s -> ", cities[bestRoute[i]].name);
    }
    printf("%s\n", cities[bestRoute[0]].name);
    printf("Best route distance: %.2f km\n", bestDistance);
}

int branch_n_bound() {
    char fileName[50];
    char startCity[50];

    printf("Enter list of cities file name: ");
    scanf("%s", fileName);
    printf("Enter starting point: ");
    scanf("%s", startCity);

    citiesFile_bnb(fileName);
    if (numCities == 0) {
        printf("No cities read. Please check the input file.\n");
        return 1;
    }
    computeDistances();

    int startIndex = findStartingCity(startCity);
    if (startIndex == -1) {
        printf("Starting city not found: %s\n", startCity);
        return 1;
    }

    int route[15];
    route[0] = startIndex;

    struct timeval start, end;
    gettimeofday(&start, NULL);
    branchbound(startIndex, 1 << startIndex, 0, route, 1);
    gettimeofday(&end, NULL);

    double elapsed = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6;

    printRoute();
    printf("Time elapsed: %.10f s\n", elapsed);

    return 0;
}
// END OF BRANCH N BOUND

//  START OF ACO
typedef struct ant {
    int order_ant[MAX_CITY]; // Urutan didatanginya kota
    int visited_ant[MAX_CITY]; // Penanda kota yang sudah dilewati
    int current_city; // Posisi semut sekarang
    int n_city; // Jumlah kota yang sudah didatangi
    double distance_ant; // Jarak total yang dilewati
} ant;

ant ants[MAX_ANT];
char cities_ant[MAX_CITY][MAX_NAME];
double latitude[MAX_CITY], longitude[MAX_CITY], distance_ant[MAX_CITY][MAX_CITY], pheromones[MAX_CITY][MAX_CITY];

// Fungsi untuk menentukan jarak antar 2 kota
double haversine_aco(double lat1, double lon1, double lat2, double lon2) {
    double dLat = (lat2 - lat1) * PI / 180.0;
    double dLon = (lon2 - lon1) * PI / 180.0;

    double a = sin(dLat / 2) * sin(dLat / 2) +
               cos(lat1 * PI / 180.0) * cos(lat2 * PI / 180.0) * 
               sin(dLon / 2) * sin(dLon / 2);
    double c = 2 * asin(sqrt(a));
    double distance_ant = 6371.0 * c;

    return distance_ant;
}

// Inisialisasi satu struct ant, berdasarkan input kota asal
int initialize_ant(ant *a, char city[]) {
    int found = -1;
    for (int i = 0; i < MAX_CITY; i++) {
        if (strcmp(city, cities_ant[i]) == 0) {
            found = i;
            break;
        }
    }

    if (found != -1) {
        for (int i = 0; i < MAX_CITY; i++) {
            a->visited_ant[i] = 0;
        }
        a->visited_ant[found] = 1;
        a->current_city = found;
        a->order_ant[0] = found;  // Initialize the starting city
        a->distance_ant = 0.0;
        a->n_city = 1;
        return 1;
    } else {
        printf("Starting point not found\n");
        return 0;
    }
}

// Fungsi untuk mengevaporasi feromon di setiap edge
void pheromone_evaporator() {
    // Evaporate pheromones first
    for (int i = 0; i < MAX_CITY; i++) {
        for (int j = 0; j < MAX_CITY; j++) {
            pheromones[i][j] *= (1 - RHO);
            if (pheromones[i][j] < 1.0) {
                pheromones[i][j] = 1.0;
            }
        }
    }
}

// Fungsi untuk menyimpan feromon di setiap edge yang dilewati semut
void pheromone_generator(ant a) {
    // Update pheromones based on the ant's tour
    for (int i = 0; i < a.n_city - 1; i++) {
        pheromones[a.order_ant[i]][a.order_ant[i + 1]] += Q / a.distance_ant;   // Jumlah feromon tergantung jarak yang dilewati
        pheromones[a.order_ant[i + 1]][a.order_ant[i]] += Q / a.distance_ant;
    }
    pheromones[a.order_ant[a.n_city - 1]][a.order_ant[0]] += Q / a.distance_ant;   
    pheromones[a.order_ant[0]][a.order_ant[a.n_city - 1]] += Q / a.distance_ant;
}

//Fungsi untuk menentukan kota selanjutnya
int next_city(ant a, int total_city) {
    double desires[total_city], sum_desire = 0.0;

    for (int i = 0; i < total_city; i++) {
        if (a.visited_ant[i] == 0) {
            desires[i] = pow(pheromones[a.current_city][i], ALPHA) * pow(1.0 / distance_ant[a.current_city][i], BETA); // Pemilihan bergantung pada jumlah feromon dan jarak setiap edge
            if (desires[i] < 1.0) {
                desires[i] = 1.0;
            }
            sum_desire += desires[i];
        } else {
            desires[i] = 0.0;
        }
    }

    // Pemilihan secara acak berbobot 
    double random = ((double)rand() / RAND_MAX) * sum_desire;
    double cumulative_probability = 0.0;

    for (int i = 0; i < total_city; i++) {
        if (desires[i] > 0) {
            cumulative_probability += desires[i];
            if (random <= cumulative_probability) {
                return i;
            }
        }
    }

    // Kasus kota tidak terpilih (seharusnya tidak mungkin terjadi)
    for (int i = 0; i < total_city; i++) {
        if (a.visited_ant[i] == 0) {
            return i;
        }
    }
    return -1;
}

// Menjalankan sebuah semut
void ant_simulator(ant *a, int total_city) {
    int city = 0;
    while (city != -1) {
        city = next_city(*a, total_city);
        if (city != -1) {
            a->visited_ant[city] = 1;  
            a->order_ant[a->n_city] = city;
            a->n_city++;
            a->distance_ant += distance_ant[a->current_city][city];
            a->current_city = city;
        }
    }
    a->distance_ant += distance_ant[a->current_city][a->order_ant[0]];
}

// Mencetak jalur yang dilewati satu semut
void print_tour(ant *a) {
    printf("Best route found:\n");
    for (int i = 0; i < a->n_city; i++) {
        printf("%s -> ", cities_ant[a->order_ant[i]]);
    }
    printf("%s\n", cities_ant[a->order_ant[0]]);
    printf("Best route distance: %lf km\n", a->distance_ant);
}

// Fungsi utama ACO
int aco() {
    char filename[20];
    char line[1000];

    printf("Enter list of cities file name: ");
    scanf("%s", filename);

    // Pembukaan file
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error opening file\n");
        return 1;
    }

    // Pembacaan file
    char *token;
    int counter = 0;
    while (fgets(line, sizeof(line), file)) {
        token = strtok(line, ",");
        strcpy(cities_ant[counter], token);
        token = strtok(NULL, ",");
        latitude[counter] = atof(token);
        token = strtok(NULL, "\n");
        longitude[counter] = atof(token);
        counter++;
    }
    fclose(file);

    for (int i = 0; i < counter; i++) {
        for (int j = 0; j < counter; j++) {
            // Perhitungan jarak setiap edge
            distance_ant[i][j] = haversine_aco(latitude[i], longitude[i], latitude[j], longitude[j]);

            //Inisialisasi nilai feromon ke nilai yang sangat kecil
            pheromones[i][j] = 0.001; 
        }
    }

    // Pengambilan input kota asal
    char initial_city[MAX_NAME];
    printf("Enter starting point: ");
    scanf("%s", initial_city);

    // Seed randomizer
    srand(time(NULL)); 

    int success;

    clock_t start_time = clock(); // Start time

    // Loop mobilisasi semut
    for (int j = 0; j < MAX_ITERATIONS; j++) {
        // Mobilisasi semut secara "paralel" 

        //inisialisasi
        for (int i = 0; i < MAX_ANT; i++) {
            success = initialize_ant(&ants[i], initial_city);
            if (!success) break;
        }

        //mobilisasi
        if (success) {
            for (int i = 0; i < MAX_ANT; i++) {
                ant_simulator(&ants[i], counter);
            }
            for (int i = 0; i < MAX_ANT; i++) {
                pheromone_generator(ants[i]);
            }
            pheromone_evaporator();
        }
        else break;
    }

    clock_t end_time = clock(); // End time
    double time_spent = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // Pencetakan hasil perjalan semut paling akhir
    if (success){
        print_tour(&ants[MAX_ANT-1]);
    
        printf("Time elapsed: %lf seconds\n", time_spent);
    }

    return 0;
}
// END OF ACO

// START OF GENETIC
typedef struct {
    char name[50];
    double latitude;
    double longitude;
} City_genetic;

int NUM_CITIES = NUM_CITIES_init;

City_genetic cities_genetic[NUM_CITIES_init];
double distances_genetic[NUM_CITIES_init][NUM_CITIES_init];
int starting_city_genetic = 0;  // User-defined starting city index

void read_cities(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    char line[100];
    int i = 0;
    while (fgets(line, sizeof(line), file) && i < NUM_CITIES) {
        char *token = strtok(line, ",");
        strcpy(cities_genetic[i].name, token);
        token = strtok(NULL, ",");
        cities_genetic[i].latitude = atof(token);
        token = strtok(NULL, ",");
        cities_genetic[i].longitude = atof(token);
        i++;
    }
    NUM_CITIES = i;
    fclose(file);
}

double calculate_distance_genetic(double lat1, double lon1, double lat2, double lon2) {
    double dlat = (lat2 - lat1) * PI / 180.0;
    double dlon = (lon2 - lon1) * PI / 180.0;
    lat1 = lat1 * PI / 180.0;
    lat2 = lat2 * PI / 180.0;

    double a = sin(dlat / 2) * sin(dlat / 2) + cos(lat1) * cos(lat2) * sin(dlon / 2) * sin(dlon / 2);
    double c = 2 * asin(sqrt(a));
    return EARTH_RADIUS * c;
}

void calculate_all_distance() {
    for (int i = 0; i < NUM_CITIES; i++) {
        for (int j = 0; j < NUM_CITIES; j++) {
            distances_genetic[i][j] = (float)calculate_distance_genetic
        (cities_genetic[i].latitude, cities_genetic[i].longitude, cities_genetic[j].latitude, cities_genetic[j].longitude);
        }
    }
}

void initialize_population(int population[POP_SIZE][NUM_CITIES]) {
    for (int i = 0; i < POP_SIZE; i++) {
        int city = 0;
        for (int j = 0; j < NUM_CITIES - 1; j++) {
            if (city == starting_city_genetic) city++;  // Skip the starting city
            population[i][j] = city++;
        }
        // Shuffle the cities_genetic to create a random tour
        for (int j = 0; j < NUM_CITIES - 2; j++) {  // Exclude the last city slot
            int r = rand() % (NUM_CITIES - 2);
            int temp = population[i][j];
            population[i][j] = population[i][r];
            population[i][r] = temp;
        }
        population[i][NUM_CITIES - 1] = starting_city_genetic;  // Ensure the starting city is at the end
    }
}

void calculate_fitness(int population[POP_SIZE][NUM_CITIES], double fitness[POP_SIZE]) {
    for (int i = 0; i < POP_SIZE; i++) {
        double total_distance = 0;
        total_distance += distances_genetic[starting_city_genetic][population[i][0]];
        for (int j = 0; j < NUM_CITIES - 2; j++) {
            total_distance += distances_genetic[population[i][j]][population[i][j + 1]];
        }
        total_distance += distances_genetic[population[i][NUM_CITIES - 2]][starting_city_genetic]; // Return to the starting city
        fitness[i] = total_distance;
    }
}

void copy_individual(int src[NUM_CITIES], int dest[NUM_CITIES]) {
    for (int i = 0; i < NUM_CITIES; i++) {
        dest[i] = src[i];
    }
}

void crossover(int parent1[NUM_CITIES], int parent2[NUM_CITIES], int child[NUM_CITIES]) {
    int start = rand() % (NUM_CITIES - 1);
    int end = rand() % (NUM_CITIES - 1);

    if (start > end) {
        int temp = start;
        start = end;
        end = temp;
    }

    for (int i = start; i <= end; i++) {
        child[i] = parent1[i];
    }

    int index = 0;
    for (int i = 0; i < NUM_CITIES - 1; i++) {
        int city = parent2[i];
        int found = 0;
        for (int j = start; j <= end; j++) {
            if (child[j] == city) {
                found = 1;
                break;
            }
        }
        if (!found) {
            while (index >= start && index <= end) {
                index++;
            }
            child[index++] = city;
        }
    }
    child[NUM_CITIES - 1] = starting_city_genetic;  // Ensure the starting city is at the end
}

void mutate(int individual[NUM_CITIES]) {
    int city1 = rand() % (NUM_CITIES - 1);
    int city2 = rand() % (NUM_CITIES - 1);
    int temp = individual[city1];
    individual[city1] = individual[city2];
    individual[city2] = temp;
}

void selection(int population[POP_SIZE][NUM_CITIES], double fitness[POP_SIZE], int new_population[POP_SIZE][NUM_CITIES]) {
    for (int i = 0; i < POP_SIZE; i++) {
        int parent1 = rand() % POP_SIZE;
        int parent2 = rand() % POP_SIZE;
        if (fitness[parent1] < fitness[parent2]) {
            copy_individual(population[parent1], new_population[i]);
        } else {
            copy_individual(population[parent2], new_population[i]);
        }

        // Crossover with another parent
        if (rand() % 2) {
            int partner = rand() % POP_SIZE;
            crossover(new_population[i], population[partner], new_population[i]);
        }

        // Mutate the child
        if ((float)rand() / RAND_MAX < MUTATION_RATE) {
            mutate(new_population[i]);
        }
    }
}

int find_best_solution(double fitness[POP_SIZE]) {
    int best_index = 0;
    for (int i = 1; i < POP_SIZE; i++) {
        if (fitness[i] < fitness[best_index]) {
            best_index = i;
        }
    }
    return best_index;
}

int find_city_index(const char *city_name) {
    for (int i = 0; i < NUM_CITIES; i++) {
        if (strcmp(cities_genetic[i].name, city_name) == 0) {
            return i; // Return the index of the city if found
        }
    }
    return -1; // Return -1 if city not found
}

int genetics() {
    srand(time(0));

    // Read cities_genetic from CSV file
    char fileName[50];
    printf("Enter list of cities file name: ");
    scanf("%s", fileName);
    read_cities(fileName);
    
    // Calculate distances_genetic between cities_genetic

    calculate_all_distance();

    // Allow user to choose the starting city by name
    char starting_city_name[50];
    printf("Enter the starting city name: ");
    scanf("%s", starting_city_name);
    starting_city_genetic = find_city_index(starting_city_name);
    clock_t begin = clock();
    
    if (starting_city_genetic == -1) {
        printf("City_genetic not found!\n");
        return 1;
    }

    int population[POP_SIZE][NUM_CITIES];
    int new_population[POP_SIZE][NUM_CITIES];
    double fitness[POP_SIZE];

    initialize_population(population);

    for (int gen = 0; gen < MAX_GEN; gen++) {
        calculate_fitness(population, fitness);
        selection(population, fitness, new_population);

        // Copy new_population to population
        for (int i = 0; i < POP_SIZE; i++) {
            copy_individual(new_population[i], population[i]);
        }
    }

    calculate_fitness(population, fitness);
    int best_index = find_best_solution(fitness);
    
    printf("Best solution found:\n");
    printf("%s-> ", cities_genetic[starting_city_genetic].name);
    for (int i = 0; i < NUM_CITIES - 1; i++) {
        if (cities_genetic[population[best_index][i]].name != "\0"){
            printf("%s-> ", cities_genetic
        [population[best_index][i]].name);
        }
    }
    clock_t end = clock();
    double timeSpent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("%s\n", cities_genetic[starting_city_genetic].name);  // Print the starting city again at the end
    printf("Total distance: %lf km\n", fitness[best_index]);
    printf("Time elapsed: %lf s\n", timeSpent);

    return 0;
}
// END OF GENETICS

int main (){
    int menu;
    printf("Selamat Datang di Program Best Algorithm!\nSilahkan pilih jenis algoritma apa yang ingin diakses.\n1. Deep First Search (DFS)\n2. Breadth First Search(BFS)\n3. Greedy Algorithm\n4. Brute Force\n5. Branch and Bound\n6. Genetic Algorithm\n7. Ant Colony Optimization\n Algoritma yang ingin diakses adalah algoritma (angka-nya saja) :");
    scanf("%d", &menu);
    if (menu == 1){
        dfs();
    }
    else if (menu == 2){
        bfs();
    }
    else if (menu == 3){
        greedy();
    }
    else if (menu == 4){
        brute_force();
    }
    else if (menu == 5){
        branch_n_bound();
    }
    else if (menu == 6){
        genetics();
    }
    else if (menu == 7){
        aco();
    }
    else{
        printf("Input salah, program akan dimatikan, silahkan start ulang ^^");
    }
}
