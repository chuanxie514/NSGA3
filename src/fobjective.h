#ifndef FOBJECTIVE_H_
#define FOBJECTIVE_H_

#include "global.h"

#define pi   3.1415926
#define SQR2 sqrt(2)

void objectives(vector<double> x_var, vector<double>& y_obj)
{
    if (!strcmp(strTestInstance, "ZDT1")) {
        y_obj[0] = x_var[0];
        double g = 0;
        for (int i = 1; i < numVariables; i++)
            g += x_var[i];
        g = 1 + 9 * g / (numVariables - 1);
        y_obj[1] = g * (1 - sqrt(x_var[0] / g));
    }
    else if (!strcmp(strTestInstance, "ZDT2")) {
        y_obj[0] = x_var[0];
        double g = 0;
        for (int i = 1; i < numVariables; i++)
            g += x_var[i];
        g = 1 + 9 * g / (numVariables - 1);
        y_obj[1] = g * (1 - pow(x_var[0] / g, 2));
    }
    else if (!strcmp(strTestInstance, "ZDT3")) {
        y_obj[0] = x_var[0];
        double g = 0;
        for (int i = 1; i < numVariables; i++)
            g += x_var[i];
        g = 1 + 9 * g / (numVariables - 1);
        y_obj[1] = g * (1 - sqrt(x_var[0] / g) - (x_var[0] / g) * sin(10 * pi * x_var[0]));
    }
    else if (!strcmp(strTestInstance, "ZDT4")) {
        y_obj[0] = x_var[0];
        double g = 0;
        for (int i = 1; i < numVariables; i++)
            g += pow(x_var[i], 2) - 10 * cos(4 * pi * x_var[i]);
        g += 1 + 10 * (numVariables - 1);
        y_obj[1] = g * (1 - sqrt(x_var[0] / g));
    }
    else if (!strcmp(strTestInstance, "ZDT6")) {
        y_obj[0] = 1 - exp(-4 * x_var[0]) * pow(sin(6 * pi * x_var[0]), 6);
        double g = 0;
        for (int i = 1; i < numVariables; i++)
            g += x_var[i];
        g = 1 + 9 * pow(g / (numVariables - 1), 0.25);
        y_obj[1] = g * (1 - pow(y_obj[0] / g, 2));
    }
    else if (!strcmp(strTestInstance, "OKA2")) {
        double x1 = x_var[0];
        double x2 = max(min(x_var[1], 5.0), -5.0);
        double x3 = max(min(x_var[2], 5.0), -5.0);

        y_obj[0] = x1;
        y_obj[1] = 1 - (1 / (4 * pi * pi)) * pow(x1 + pi, 2) + cbrt(fabs(x2 - 5 * cos(x1))) + cbrt(fabs(x3 - 5 * sin(x1)));
    }
    else if (!strcmp(strTestInstance, "DTLZ1")) {
        int k = numVariables - numObjectives + 1;
        double g = 0;
        for (int i = numVariables - k; i < numVariables; i++)
            g += pow(x_var[i] - 0.5, 2) - cos(20 * pi * (x_var[i] - 0.5));
        g = 100 * (k + g);

        for (int i = 0; i < numObjectives; i++) {
            double f = 0.5 * (1 + g);
            for (int j = 0; j < numObjectives - i - 1; j++)
                f *= x_var[j];
            if (i != 0)
                f *= 1 - x_var[numObjectives - i - 1];
            y_obj[i] = f;
        }
    }
    else if (!strcmp(strTestInstance, "DTLZ2")) {
        int k = numVariables - numObjectives + 1;
        double g = 0;
        for (int i = numVariables - k; i < numVariables; i++)
            g += pow(x_var[i] - 0.5, 2);

        for (int i = 0; i < numObjectives; i++) {
            double f = (1 + g);
            for (int j = 0; j < numObjectives - i - 1; j++)
                f *= cos(x_var[j] * 0.5 * pi);
            if (i != 0)
                f *= sin(x_var[numObjectives - i - 1] * 0.5 * pi);
            y_obj[i] = f;
        }
    }
    else if (!strcmp(strTestInstance, "DTLZ3")) {
        int k = numVariables - numObjectives + 1;
        double g = 0;
        for (int i = numVariables - k; i < numVariables; i++)
            g += pow(x_var[i] - 0.5, 2) - cos(20 * pi * (x_var[i] - 0.5));
        g = 100 * (k + g);

        for (int i = 0; i < numObjectives; i++) {
            double f = (1 + g);
            for (int j = 0; j < numObjectives - i - 1; j++)
                f *= cos(pow(x_var[j], 100) * 0.5 * pi);
            if (i != 0)
                f *= sin(pow(x_var[numObjectives - i - 1], 100) * 0.5 * pi);
            y_obj[i] = f;
        }
    }
    else if (!strcmp(strTestInstance, "DTLZ4")) {
        int k = numVariables - numObjectives + 1;
        double g = 0;
        for (int i = numVariables - k; i < numVariables; i++)
            g += pow(x_var[i] - 0.5, 2);

        for (int i = 0; i < numObjectives; i++) {
            double f = 1 + g;
            for (int j = 0; j < numObjectives - i - 1; j++)
                f *= cos(pow(x_var[j], 100) * 0.5 * pi);
            if (i != 0)
                f *= sin(pow(x_var[numObjectives - i - 1], 100) * 0.5 * pi);
            y_obj[i] = f;
        }
    }
    else if (!strcmp(strTestInstance, "DTLZ5")) {
        int k = numVariables - numObjectives + 1;
        double g = 0;
        for (int i = numVariables - k; i < numVariables; i++)
            g += pow(x_var[i] - 0.5, 2);

        vector<double> theta(numObjectives - 1);
        for (int i = 0; i < numObjectives - 1; i++)
            theta[i] = pi / 2 * (1 + 2 * g * x_var[i]);

        for (int i = 0; i < numObjectives; i++) {
            double f = 1 + g;
            if (i != 0)
                f *= cos(theta[i - 1]);
            if (i != numObjectives - 1)
                f *= sin(theta[i]);
            y_obj[i] = f;
        }
    }
    else if (!strcmp(strTestInstance, "DTLZ6")) {
        int k = numVariables - numObjectives + 1;
        double g = 0;
        for (int i = numVariables - k; i < numVariables; i++)
            g += pow(x_var[i], 0.1);

        vector<double> theta(numObjectives - 1);
        for (int i = 0; i < numObjectives - 1; i++)
            theta[i] = pi / 2 * (1 + 2 * g * x_var[i]);

        for (int i = 0; i < numObjectives; i++) {
            double f = 1 + g;
            if (i != 0)
                f *= cos(theta[i - 1]);
            if (i != numObjectives - 1)
                f *= sin(theta[i]);
            y_obj[i] = f;
        }
    }
}

#endif /* FOBJECTIVE_H_ */
