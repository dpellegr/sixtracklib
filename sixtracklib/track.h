typedef struct {
    double length;
} Drift ;

typedef struct {
      double length;
} DriftExact ;

typedef struct {
      long int order;
      double l ;
      double hxl;
      double hyl;
      double bal[1];
} Multipole;

typedef struct {
      double volt;
      double freq;
      double lag;
} Cavity;

typedef struct {
      double cz;
      double sz;
      double dx;
      double dy;
} Align;
