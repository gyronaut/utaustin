
void rootTest(){
    //declaring basic variables (int), assigning them values
    int a = 5;
    int b = 7;

    //declaring a new value, assinging it a value based on other variables
    int c = a+b;

    /* Printing out a statement to the command line:
     *
     * the function printf() takes a string as a parameter and prints it on the command line.  In the second case, I used "%i"
     * as a placeholder for an int, and then told it to fill that placeholder with the value stored in the variable
     * c.
     *
     * The last print statement shows a more complicated version of this, using the values for a,b, and (a/b). Note that since
     * a and b are integers, if i hadn't included the "double", the value of a/b would have been 5/7 rounded down to an integer (0).
     *
     * Note: "\n" is the new line character
     */

    printf("Printing to the terminal.\n");

    printf("a+b = %i \n", c);

    printf("a=%i, b=%i, a/b = %f \n", a, b, double(a)/double(b));
    
    
    /* Declaring a new TH1F histogram:
     * 
     * this declares a new histogram histo with 100 bins between x=a and x=b
     */

    TH1F histo("histo", "Histogram", 100, a,b);

    /* Declaring a new TH1F histogram via pointer:
     *
     * First thing to note is the * on TH1F*.  This means we're defining a "pointer" to a histogram, and not a histogram itself.
     * The declaration is saying "this pointer named 'histoPointer' is the memory address of a new histogram with 20 bins between x=a and x=b".
     *
     * Whole books can (and have) been written about pointers, but the important parts are:
     *
     *** "histoPointer" isn't a histogram, it's a memory address
     *** "*histoPointer" is the histogram itself (the '*' says "the object this pointer is pointing to")
     *** "histo" is a histogram
     *** "&histo" is the memory address where 'histo' is stored (the '&' says "the memory address of this object")
     */
    TH1F* histoPointer = new TH1F("ctest","ctest", 20, a, b);


    //calling functions of a histogram
    //Note that functions of an object are called with a "."
    histo.FillRandom("gaus", 10000);

    //Declaring a new canvas (pointer):
    TCanvas* c1 = new TCanvas("canvas", "canvas", 800, 1000);
    c1->cd();

    //calling functions of a pointer to a histogram:
    //Note that functions of a POINTER to an object are called with a "->"
    histoPointer->FillRandom("gaus", 500);
    histoPointer->Draw();
}
