#ifndef ERCounter_hpp
#define ERCounter_hpp

class ERCounter
{
    protected:
        int errors = 0;                                                         // Number of symbol errors.
        int corrects = 0;                                                       // Number of correct bits.

    public:
        ERCounter();                                                            // Constructor.
        virtual ~ERCounter();                                                   // Deconstructor.
        void clear();                                                           // Clear counts.
        int get_errors() const;                                                 // Get the number of symbol errors.
        int get_corrects() const;                                               // Get the number of correct bits.
        double get_errorrate() const;                                           // Get symbol error rate.
};

#endif
