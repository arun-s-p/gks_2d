#include <string>

// Assuming these variables are global or defined in another header file
extern int IREAD, IOUT, ISIM, IFLO;
extern std::string FNAME0, FNAME1, FNAME2;

// Equivalent structure to encapsulate input/output definitions
struct IODefinitions {
    int IREAD;
    int IOUT, ISIM, IFLO;
    std::string FNAME0, FNAME1, FNAME2;
};

// Function to initialize input/output definitions
IODefinitions InitializeIODefinitions() {
    IODefinitions ioDefs;
    ioDefs.IREAD = 0;
    ioDefs.IOUT = 0;
    ioDefs.ISIM = 0;
    ioDefs.IFLO = 0;
    ioDefs.FNAME0 = "";
    ioDefs.FNAME1 = "";
    ioDefs.FNAME2 = "";

    return ioDefs;
}

