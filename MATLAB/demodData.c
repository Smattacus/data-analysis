#include "mex.h"
#include "math.h"

#define TRUE 1
#define FALSE 0
#define MAX_COL_LENGTH 10

void countTopsBots(double *, int , int *, int *);
void demodData(double *, double *, int, double *, double *) ;

void mexFunction(int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{
	/*Data inputs: prhs[0] = Data
	prhs[1] = Square wave.

	Data Outputs:
	plhs[0] = Tops
	prhs[1] = Bottoms

EG: [T, B] = demodData(data, square_wave)

	*/

	size_t mrows_d, ncols_d, mrows_sq, ncols_sq;
	size_t sq_length, d_length; int tops, bots;
	double *sq; double *data;
	double *top_pr, *bot_pr;
	mxArray *out_top, *out_bot;

/* Check for proper number of arguments */
	if (nrhs != 2) {
		mexErrMsgIdAndTxt("MATLAB:demodData:invalidNumInputs",
			"Two inputs required.");
	} else if(nlhs != 2) {
		mexErrMsgIdAndTxt("MATLAB:demodData:invalidNumInputs",
			"Two outputs required.");
	}

	/*The inputs should be one dimensional arrays. */
	mrows_d = mxGetM(prhs[0]);
	ncols_d = mxGetN(prhs[0]);
	mrows_sq = mxGetM(prhs[1]);
	ncols_sq = mxGetN(prhs[1]);
	//Check that the input arrays are 1xN or Mx1.
	if ((mrows_d != 1 && ncols_d != 1))mexErrMsgIdAndTxt("MATLAB:demodData:invalidNumInputs",
		"Data Input must be 1-d array.");
	if (mrows_d == 1 && ncols_d == 1) mexErrMsgIdAndTxt("MATLAB:demodData:invalidNumInputs",
		"Data array should not be a scalar.");
	if ((mrows_sq != 1 && ncols_sq != 1))mexErrMsgIdAndTxt("MATLAB:demodData:invalidNumInputs",
		"Square Wave input must be 1-d array.");
	if (mrows_sq == 1 && ncols_sq == 1) mexErrMsgIdAndTxt("MATLAB:demodData:invalidNumInputs",
		"Square wave array should not be a scalar.");

	//Load the non-singleton dimension.
	if (mrows_sq == 1) sq_length = ncols_sq; else sq_length = mrows_sq;
	if (mrows_d == 1) d_length = ncols_d; else d_length = mrows_d;

	//Get the square data array pointer.
	sq = mxGetPr(prhs[1]);
	data = mxGetPr(prhs[0]);
	//Count the tops and bottoms of the square wave.
	countTopsBots(sq, d_length, &tops, &bots);

	//Create the two output arrays
	out_top = mxCreateDoubleMatrix(MAX_COL_LENGTH, tops, mxREAL);
	out_bot = mxCreateDoubleMatrix(MAX_COL_LENGTH, bots, mxREAL);
	top_pr = mxGetPr(out_top);
	bot_pr = mxGetPr(out_bot);

	//Populate the two output arrays.
	demodData(data, sq, d_length, top_pr, bot_pr);
	plhs[0] = out_top;
	plhs[1] = out_bot;

}


void countTopsBots(double *sq, int length, int *tops, int *bots)
{
	/**Helper function to count the number of tops and bottoms. */
	int i; int top_bool;
	int temp_top = 0; int temp_bot = 0;
	if (sq[0] == 0) { top_bool = FALSE; temp_bot++;} else { top_bool = TRUE; temp_top++;}
	for (i=1; i< length;i++)
	{
		//Check if the square wave has made a jump - continue if it hasn't.
		if ((top_bool == TRUE && sq[i] == 2) || (top_bool == FALSE && sq[i] == 0)) continue;
		//Switch top_bool and increment the counter.
		if (top_bool == TRUE) temp_top++; else temp_bot++;
		top_bool = !top_bool;
	}
	if (temp_top > temp_bot) temp_bot = temp_top; else temp_top = temp_bot;
	*tops = temp_top; *bots = temp_bot;
}

void demodData(double *data, double *square, int length, double *top, double *bot) {

	int i; int top_bool; int bot_count = 0; int top_count=0; int rc =0;

	//Initialize the booleans before jumping into the for loop
	if (square[0] == 0) {
		top_bool = FALSE;// bot[0] = data[0]; rc++;
	} else { 
		top_bool = TRUE;// top[0] = data[0]; rc++;
	} 

	//Step through the entire data array, checking jumps on square
	//to determine time to switch array being loaded.
	for(i=0; i < length; i++)
	{
		if (top_bool == TRUE && square[i] == 2) {
			top[rc++ + top_count * MAX_COL_LENGTH] = data[i];
			continue;
		} else if (top_bool == FALSE && square[i] == 0) {
			bot[rc++ + bot_count * MAX_COL_LENGTH] = data[i];
			continue;
		} else if (top_bool == FALSE && square[i] != 0) {
			top_bool = !top_bool;
			//Reset row counter to 0, increment column count.
			rc = 0; bot_count++;
			top[rc++ + top_count * MAX_COL_LENGTH] = data[i];
			continue;
		} else if (top_bool == TRUE && square[i] != 2) {
			top_bool = !top_bool; rc =0; top_count++;
			bot[rc++ + bot_count * MAX_COL_LENGTH] = data[i];
			continue;
		}
	}
}