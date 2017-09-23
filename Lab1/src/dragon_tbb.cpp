/*
 * dragon_tbb.c
 *
 *  Created on: 2011-08-17
 *      Author: Francis Giraldeau <francis.giraldeau@gmail.com>
 */

#include <iostream>

extern "C" {
#include "dragon.h"
#include "color.h"
#include "utils.h"
}
#include "dragon_tbb.h"
#include "tbb/tbb.h"
#include "TidMap.h"

using namespace std;
using namespace tbb;

static uint64_t* thread_limits;
static int nb_threads;
static uint64_t range_size;

class DragonLimits {
private:
	piece_t master;

public:
	const piece_t& getMaster() 
	{
		return master;
	}

	void operator() (const blocked_range<uint64_t> &r)
	{
		uint64_t start = r.begin();
		uint64_t end   = r.end();
		piece_limit(start, end, &master);
	}

	void join(DragonLimits &l)
	{
		piece_merge(&master, l.master);
	}

	DragonLimits(DragonLimits &x, split)
	{
		piece_init(&master);
	}

	DragonLimits()
	{
		piece_init(&master);
	}
};

class DragonDraw {
private:
	struct draw_data data;
	TidMap *tidMap;
public:
	void operator() (const blocked_range<uint64_t> &r) const
	{
		/* Go through the list `data.tid`,
		   Find the first index at which the value is 0,
		   Use that value as thread index */
		// int id = 0;
		// while (data.tid[id] != 0) id++;
		// data.tid[id] = 1;	// set this thread to BUSY.
		//int id = tidMap->getIdFromTid(gettid());
		tidMap->getIdFromTid(gettid());
		int color_id = r.begin()*range_size/nb_threads;
		//while(r.end()>=thread_limits[color_id+1]) {
			//color_id++;
		//}

		
		dragon_draw_raw(r.begin(), r.end(), data.dragon, data.dragon_width, data.dragon_height, data.limits, color_id);
		// data.tid[id] = 0;	// set this thread to FREED.
	}

	DragonDraw(struct draw_data &data, TidMap* tidMap):data(data), tidMap(tidMap) {}
};

class DragonRender {
private:
	struct draw_data data;
public:
	void operator() (const blocked_range<int> &r) const
	{
		/**/
		scale_dragon(r.begin(), r.end(), data.image, data.image_width, data.image_height, data.dragon, data.dragon_width, data.dragon_height, data.palette);
	}

	DragonRender(struct draw_data &data):data(data) {}
};

class DragonClear {
private:
	struct draw_data data;

public:

	void operator() (const blocked_range<int> &r) const
	{
		init_canvas(r.begin(), r.end(), data.dragon, -1);
	}

	DragonClear(struct draw_data &data):data(data) {}
};

int dragon_draw_tbb(char **canvas, struct rgb *image, int width, int height, uint64_t size, int nb_thread)
{
	//TODO("dragon_draw_tbb");
	TidMap tidMap(nb_thread);
	struct draw_data data;
	limits_t limits;
	char *dragon = NULL;
	int dragon_width;
	int dragon_height;
	int dragon_surface;
	int scale_x;
	int scale_y;
	int scale;
	int deltaJ;
	int deltaI;

	thread_limits = (uint64_t*)malloc(sizeof(uint64_t) * (nb_thread+1));
	nb_threads = nb_thread;
	for(int j=0; j < nb_thread ; ++j) {
		thread_limits[j] = j*size/nb_thread;
	}
	thread_limits[nb_thread] = size;
	range_size = size;

	struct palette *palette = init_palette(nb_thread);
	if (palette == NULL)
		return -1;

	/* 1. Calculer les limites du dragon */
	dragon_limits_tbb(&limits, size, nb_thread);

	task_scheduler_init init(nb_thread);

	dragon_width = limits.maximums.x - limits.minimums.x;
	dragon_height = limits.maximums.y - limits.minimums.y;
	dragon_surface = dragon_width * dragon_height;
	scale_x = dragon_width / width + 1;
	scale_y = dragon_height / height + 1;
	scale = (scale_x > scale_y ? scale_x : scale_y);
	deltaJ = (scale * width - dragon_width) / 2;
	deltaI = (scale * height - dragon_height) / 2;

	dragon = (char *) malloc(dragon_surface);
	if (dragon == NULL) {
		free_palette(palette);
		return -1;
	}

	data.nb_thread = nb_thread;
	data.dragon = dragon;
	data.image = image;
	data.size = size;
	data.image_height = height;
	data.image_width = width;
	data.dragon_width = dragon_width;
	data.dragon_height = dragon_height;
	data.limits = limits;
	data.scale = scale;
	data.deltaI = deltaI;
	data.deltaJ = deltaJ;
	data.palette = palette;
	data.tid = (int *) calloc(nb_thread, sizeof(int));

	/* 2. Initialiser la surface : DragonClear */
	DragonClear clear{data};
	parallel_for( blocked_range<int>(0, dragon_surface), clear );

	/* 3. Dessiner le dragon : DragonDraw */
	DragonDraw draw{data, &tidMap};
	parallel_for( blocked_range<uint64_t>(0, size), draw );

	/* 4. Effectuer le rendu final */
	DragonRender render{data};
	parallel_for( blocked_range<int>(0, height), render );
	
	init.terminate();

	tidMap.dump();
	free_palette(palette);
	FREE(data.tid);
	*canvas = dragon;
	//*canvas = NULL;
	return 0;
}

/*
 * Calcule les limites en terme de largeur et de hauteur de
 * la forme du dragon. Requis pour allouer la matrice de dessin.
 */
int dragon_limits_tbb(limits_t *limits, uint64_t size, int nb_thread)
{
	//TODO("dragon_limits_tbb");
	DragonLimits lim;

	/* Create the scheduler with `nb_thread` number of threads */
	task_scheduler_init init(nb_thread);

	/* Parallelize using tbb */
	parallel_reduce( blocked_range<uint64_t>(0, size), lim );


	/* Rerive limit values */
	*limits = lim.getMaster().limits;

	return 0;
}
