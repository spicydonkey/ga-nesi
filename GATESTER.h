// A library of Mathematical functions to test a genetic algorithm
// <cmath> library must be imported

float schwefel (float *argv, int dim)
{
	float value = 418.9829f*(float)dim;
	for (int i=0; i<dim; i++)
	{
		value -= argv[i]*sin(sqrt(fabs(argv[i])));
	}
	return value;
}

double schwefel (std::vector<double> vals)
{
	int dim = vals.size();
	double value = 418.9829*(double)dim;
	
	double tmp;	
	for (int i=0; i<dim; i++)
	{
		tmp = vals[i] - 500.0;	
		value -= tmp*sin(sqrt(fabs(tmp)));
	}
	return value;
}
