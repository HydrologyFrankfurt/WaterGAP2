#include <stdio.h>


/* interpolation of temperature values */
/* with cubic splines according to     */
/* Numerical Recipes in C,             */
/* 2nd edition, chapter 3.3            */
float interp_day_temp(int n, short day, short int G_temperature[][12])
{
	const short int mid_of_month[18] = { -76, -45, -15, 16, 46, 75, 105, 136, 166,
		197, 228, 258, 289, 319, 350, 381, 411, 440
	};
	short int monthly_temp[18];
	float m1[18], m2[18], mv1[18], u[18], y[18];
	float sig, p, q;
	short int i, j;

	/* calculate spline function */

	for (i = 0; i <= 11; i++)
		monthly_temp[i + 3] = G_temperature[n][i];
	monthly_temp[0] = G_temperature[n][9];
	monthly_temp[1] = G_temperature[n][10];
	monthly_temp[2] = G_temperature[n][11];
	monthly_temp[15] = G_temperature[n][0];
	monthly_temp[16] = G_temperature[n][1];
	monthly_temp[17] = G_temperature[n][2];

	m1[0] = mid_of_month[1] - mid_of_month[0];
	mv1[0] = monthly_temp[1] - monthly_temp[0];
	for (i = 1; i <= 16; i++) {
		m1[i] = mid_of_month[i + 1] - mid_of_month[i];
		m2[i] = m1[i] + m1[i - 1];
		mv1[i] = monthly_temp[i + 1] - monthly_temp[i];
	}
	u[0] = (3. / m1[0]) * mv1[0] / m1[0];
	y[0] = -.5;
	for (i = 1; i <= 16; i++) {
		sig = m1[i - 1] / m2[i];
		p = sig * y[i - 1] + 2.0;
		y[i] = (sig - 1.0) / p;
		u[i] = (6. * (mv1[i] / m1[i] - mv1[i - 1] / m1[i - 1]) / m2[i] - sig * u[i - 1]) / p;
	}
	u[17] = (-3. / m1[16]) * mv1[16] / m1[16];
	y[17] = (u[17] - 0.5 * u[16]) / (0.5 * y[16] + 1.);
	for (i = 16; i >= 0; i--)
		y[i] = y[i] * y[i + 1] + u[i];


	/* interpolate values from january, the 1st to january, the 15th */
	if (day <= 15) {
		i = day;
		p = (mid_of_month[3] - i) / m1[2];
		q = 1 - p;
		return p * monthly_temp[2] + q * monthly_temp[3] +
			((p * p * p - p) * y[2] + (q * q * q - q) * y[3]) * (m1[2] * m1[2]) / 6.;
	}
	/* interpolate values from mid of december to the end of the year */
	if (day >= 351) {
		i = day - 350;
		q = i / m1[14];
		p = 1 - q;
		return p * monthly_temp[14] + q * monthly_temp[15] +
			((p * p * p - p) * y[14] + (q * q * q - q) * y[15]) * (m1[14] * m1[14]) / 6.;
	}

	/* interpolate values from january, the 16th to mid of december */
	for (j = 3; j <= 13; j++) {
		if (day == mid_of_month[j])
			return monthly_temp[j];
		if (day < mid_of_month[j + 1]) {
			i = day - mid_of_month[j];
			q = i / m1[j];
			p = 1 - q;
			return p * monthly_temp[j] + q * monthly_temp[j + 1] +
				((p * p * p - p) * y[j] + (q * q * q - q) * y[j + 1]) * (m1[j] * m1[j]) / 6.;
		}
	}
	return -9999;	/* never reached */
}



/*                                                         */
/* interpolate cloudiness                                  */
/* same routine as for temperature but with the constraint */
/* that daily values lies between monthly values           */
/*                                                         */
template < class T > float interp_day_sunshine(int n, short day, T G_cloud[][12])
{
	const short int mid_of_month[18] = { -76, -45, -15, 16, 46, 75, 105, 136, 166,
		197, 228, 258, 289, 319, 350, 381, 411, 440
	};
	T monthly_cloud[18];
	float m1[18], m2[18], mv1[18], u[18], y[18];
	float sig, p, q;
	short int i, j;
	const short int n_days_per_year = 365;
	T daily_cloud[n_days_per_year];

	/* calculate spline function */

	for (i = 0; i <= 11; i++)
		monthly_cloud[i + 3] = G_cloud[n][i];
	monthly_cloud[0] = G_cloud[n][9];
	monthly_cloud[1] = G_cloud[n][10];
	monthly_cloud[2] = G_cloud[n][11];
	monthly_cloud[15] = G_cloud[n][0];
	monthly_cloud[16] = G_cloud[n][1];
	monthly_cloud[17] = G_cloud[n][2];

	m1[0] = mid_of_month[1] - mid_of_month[0];
	mv1[0] = monthly_cloud[1] - monthly_cloud[0];
	for (i = 1; i <= 16; i++) {
		m1[i] = mid_of_month[i + 1] - mid_of_month[i];
		m2[i] = m1[i] + m1[i - 1];
		mv1[i] = monthly_cloud[i + 1] - monthly_cloud[i];
	}
	u[0] = (3. / m1[0]) * mv1[0] / m1[0];
	y[0] = -.5;
	for (i = 1; i <= 16; i++) {
		sig = m1[i - 1] / m2[i];
		p = sig * y[i - 1] + 2.0;
		y[i] = (sig - 1.0) / p;
		u[i] = (6. * (mv1[i] / m1[i] - mv1[i - 1] / m1[i - 1]) / m2[i] - sig * u[i - 1]) / p;
	}
	u[17] = (-3. / m1[16]) * mv1[16] / m1[16];
	y[17] = (u[17] - 0.5 * u[16]) / (0.5 * y[16] + 1.);
	for (i = 16; i >= 0; i--)
		y[i] = y[i] * y[i + 1] + u[i];

	/* interpolate values from january, the 1st to january, the 15th */
	if (day <= 15) {
		i = day;
		p = (mid_of_month[3] - i) / m1[2];
		q = 1 - p;
		daily_cloud = p * monthly_cloud[2] + q * monthly_cloud[3]
			+ ((p * p * p - p) * y[2] + (q * q * q - q) * y[3]) * (m1[2] * m1[2]) / 6.;
		if ((daily_cloud > monthly_cloud[2])
			&& (daily_cloud[day] > monthly_cloud[3])) {
			if (monthly_cloud[2] > monthly_cloud[3])
				return monthly_cloud[2];
			else
				return monthly_cloud[3];
		}
		if ((daily_cloud < monthly_cloud[2])
			&& (daily_cloud < monthly_cloud[3])) {
			if (monthly_cloud[2] < monthly_cloud[3])
				return monthly_cloud[2];
			else
				return monthly_cloud[3];
		}
		return daily_cloud;
	}

	/* interpolate values from january, the 17th to mid of december */
	if (day >= 351) {
		i = day - 350;
		q = i / m1[j];
		p = 1 - q;
		daily_cloud = p * monthly_cloud[j] + q * monthly_cloud[j + 1]
			+ ((p * p * p - p) * y[j] + (q * q * q - q) * y[j + 1]) * (m1[j] * m1[j]) / 6.;
		if ((daily_cloud > monthly_cloud[j])
			&& (daily_cloud > monthly_cloud[j + 1])) {
			if (monthly_cloud[j] > monthly_cloud[j + 1])
				return monthly_cloud[j];
			else
				return monthly_cloud[j + 1];
		}
		if ((daily_cloud < monthly_cloud[j])
			&& (daily_cloud < monthly_cloud[j + 1])) {
			if (monthly_cloud[j] < monthly_cloud[j + 1])
				return monthly_cloud[j];
			else
				return monthly_cloud[j + 1];
		}
		return daily_cloud;
	}


	/* interpolate values from mid of december to the end of the year */
	for (j = 3; j <= 13; j++) {
		if (day == mid_of_month[j])
			//return monthly_temp[j];    // stopped by mehedi due to compilation error (monthly_temp was not declared in this scope)
			return monthly_cloud[j];     // added by mehedi as a test case
		if (day < mid_of_month[j + 1]) {
			i = day - mid_of_month[j];
			q = i / m1[14];
			p = 1 - q;
			daily_cloud = p * monthly_cloud[14] + q * monthly_cloud[15]
				+ ((p * p * p - p) * y[14] + (q * q * q - q) * y[15]) * (m1[14] * m1[14]) / 6.;
			if ((daily_cloud > monthly_cloud[14])
				&& (daily_cloud > monthly_cloud[15])) {
				if (monthly_cloud[14] > monthly_cloud[15])
					return monthly_cloud[14];
				else
					return monthly_cloud[15];
			}
			if ((daily_cloud < monthly_cloud[14])
				&& (daily_cloud < monthly_cloud[15])) {
				if (monthly_cloud[14] < monthly_cloud[15])
					return monthly_cloud[14];
				else
					return monthly_cloud[15];
			}
			return daily_cloud;
		}
	}
	return -9999;	/* never reached */
}

/*                                       */
/* linear interpolation of precipitation */
/*                                       */
float interp_day_prec(int n, short day_in_month, short month, short int G_precipitation[][12])
{
	const short int n_days_per_month[12]
	= { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	const float mid_of_month_1[12]
	= { 16, 14.5, 16, 15.5, 16, 15.5, 16, 16, 15.5, 16, 15.5, 16 };

	float daily_prec;

	if ((0 == month) && (day_in_month < mid_of_month_1[0])) {
		daily_prec = (G_precipitation[n][11] / (float) n_days_per_month[11])
			+ (n_days_per_month[11] - mid_of_month_1[11] + day_in_month) *
			(((G_precipitation[n][0] / (float) n_days_per_month[0])
			  - (G_precipitation[n][11] / (float) n_days_per_month[11]))
			 / (float) (n_days_per_month[11]));
	} else {
		if ((11 == month) && (day_in_month >= mid_of_month_1[11])) {
			daily_prec = (G_precipitation[n][11] / (float) n_days_per_month[11])
				+ (day_in_month - mid_of_month_1[11]) *
				(((G_precipitation[n][0] / (float) n_days_per_month[0])
				  - (G_precipitation[n][11] / (float) n_days_per_month[11]))
				 / (float) (n_days_per_month[11]));
		} else {
			if (day_in_month >= mid_of_month_1[month]) {
				daily_prec = (G_precipitation[n][month] / (float) n_days_per_month[month])
					+ (day_in_month - mid_of_month_1[month]) *
					(((G_precipitation[n][month + 1] / (float) n_days_per_month[month + 1])
					  - (G_precipitation[n][month] / (float) n_days_per_month[month]))
					 / (float) (n_days_per_month[month]));
			} else {
				daily_prec = (G_precipitation[n][month - 1] / (float) n_days_per_month[month - 1])
					+ (n_days_per_month[month - 1]
					   - mid_of_month_1[month - 1] + day_in_month) *
					(((G_precipitation[n][month] / (float) n_days_per_month[month])
					  - (G_precipitation[n][month - 1]
						 / (float) n_days_per_month[month - 1]))
					 / (float) (n_days_per_month[month - 1]));
			}
		}
	}

	/* sometimes (very rare cases) the interpolation algorithm */
	/* leads to values below zero                              */
	/* this should be avoided                                  */
	if (daily_prec < 0)
		return 0;
	else
		return daily_prec;
}
