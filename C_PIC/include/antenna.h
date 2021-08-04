/*
    Definitions for the PFRC-1 Antenna
    Each antenna composed of 28 segments, each segment is defined over a certain region of space
*/
typedef struct
{
    double xlim[2];
    double ylim[2];
    double zlim[2];
}Segment;


typedef struct
{
    Segment* segments;
    int segment_count;
}Antenna;

extern Antenna* get_segments(double* segx_min, double* segx_max,
                             double* segy_min, double* segy_max,
                             double* segz_min, double* segz_max,
                             int num_segments);
