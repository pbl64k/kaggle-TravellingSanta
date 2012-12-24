<?php

	error_reporting(E_ALL | E_STRICT);

	ini_set('memory_limit', '6144M');

	gc_enable();

	class KdTreeNode
	{
		private $pt;
		private $data;
		private $lc = NULL;
		private $rc = NULL;

		final static public function make(array $pt, $data)
		{
			return new self($pt, $data);
		}

		final public function getPt()
		{
			return $this->pt;
		}

		final public function getData()
		{
			return $this->data;
		}

		final public function setLc($lc)
		{
			$this->lc = $lc;

			return $this;
		}

		final public function getLc()
		{
			return $this->lc;
		}

		final public function setRc($rc)
		{
			$this->rc = $rc;

			return $this;
		}

		final public function getRc()
		{
			return $this->rc;
		}

		final private function __construct(array $pt, $data)
		{
			$this->pt = $pt;
			$this->data = $data;
			$this->lc = NULL;
			$this->rc = NULL;
		}
	}

	class KdTree
	{
		private $dim;
		private $root = NULL;

		final static public function make($dim)
		{
			return new self($dim);
		}

		final public function add(array $pt, $data)
		{
			$this->root = $this->add0($pt, $data, $this->root, 0);

			return $this;
		}

		final public function findKnn(array $pt, $k)
		{
			return $this->findKnn0($pt, $k, $this->root, 0, array());
		}

		final private function __construct($dim)
		{
			$this->dim = intval($dim);
			$this->root = NULL;
		}

		final private function add0(array $pt, $data, $node, $ix)
		{
			if (is_null($node))
			{
				return KdTreeNode::make($pt, $data);
			}

			$pt0 = $node->getPt();

			if ($pt[$ix] < $pt0[$ix])
			{
				$node->setLc($this->add0($pt, $data, $node->getLc(), ($ix + 1) % $this->dim));
			}
			else
			{
				$node->setRc($this->add0($pt, $data, $node->getRc(), ($ix + 1) % $this->dim));
			}

			return $node;
		}

		final private function findKnn0(array $pt, $k, $node, $ix, array $res)
		{
			if (is_null($node))
			{
				return $res;
			}

			$res = $this->injectNn($pt, $node->getPt(), $k, $node->getData(), $res);

			$pt0 = $node->getPt();

			if ($pt[$ix] < $pt0[$ix])
			{
				$res = $this->findKnn0($pt, $k, $node->getLc(), ($ix + 1) % $this->dim, $res);

				if (count($res))
				{
					list($d, $p, $dt) = $res[count($res) - 1];
				}

				if ((count($res) < $k) || ($d > abs($pt0[$ix] - $pt[$ix])))
				{
					$res = $this->findKnn0($pt, $k, $node->getRc(), ($ix + 1) % $this->dim, $res);
				}
			}
			else
			{
				$res = $this->findKnn0($pt, $k, $node->getRc(), ($ix + 1) % $this->dim, $res);

				if (count($res))
				{
					list($d, $p, $dt) = $res[count($res) - 1];
				}

				if ((count($res) < $k) || ($d > abs($pt0[$ix] - $pt[$ix])))
				{
					$res = $this->findKnn0($pt, $k, $node->getLc(), ($ix + 1) % $this->dim, $res);
				}
			}

			return $res;
		}

		final private function injectNn(array $pt0, array $pt, $k, $data, array $res0)
		{
			global $city0, $cities, $edges, $citycheck;

			if (($citycheck && $cities[$data][2]) || isset($edges[$city0.'-'.$data]))
			{
				return $res0;
			}

			$dist = $this->calcDist($pt0, $pt);

			if (! count($res0))
			{
				return array(array($dist, $pt, $data));
			}

			$res = array();

			$injected = FALSE;

			$rsc = count($res0);
			$rc = min($k, $rsc + 1);

			for ($i = 0; $i != $rc; ++$i)
			{
				if (! $injected)
				{
					if ($i === $rsc)
					{
						$res[] = array($dist, $pt, $data);

						$injected = TRUE;

						continue;
					}

					list($d, $p, $dt) = $res0[$i - ($injected ? 1 : 0)];

					if ($dist < $d)
					{
						$res[] = array($dist, $pt, $data);

						$injected = TRUE;

						continue;
					}
				}

				$res[] = $res0[$i - ($injected ? 1 : 0)];
			}

			return $res;
		}

		final private function calcDist(array $pt0, array $pt)
		{
			return sqrt(array_sum(array_map(
					function($a, $b)
					{
						$c = $a - $b;
						return $c * $c;
					}, $pt0, $pt)));
		}
	}

	$l = explode("\n", file_get_contents(__DIR__.'/santa_cities.csv'));

	//$bucketsize = 500;
	//$bucketsize = 1000;
	//$bucketsize = 2000;
	$bucketsize = 20000;
	$bnum = 20000 / $bucketsize;

	//$runs = 3;
	//$runs = 5;
	//$runs = 40;
	$runs = 1;

	$kc = 16;
	$k0 = 2;
	
	$opt = 150000;
	//$opt = 6000;

	$showPaths = TRUE;

	$read = TRUE;

	$citycheck = TRUE;

	$fn = __DIR__.'/sss-dat-'.$bucketsize.'.ser';

	if ($read)
	{
		list($cities, $buckets, $kdts) = unserialize(file_get_contents($fn));
	}
	else
	{
		$buckets = array();
		$kdts = array();
	
		$cities = array();
	
		foreach ($l as $n => $line)
		{
			if (! $n)
			{
				continue;
			}
	
			$x = array_map('intval', explode(',', $line));
	
			if (count($x) < 3)
			{
				break;
			}
	
			print('Importing line '.$n.'...'."\n");
	
			$xbucket = intval(floor($x[1] / $bucketsize));
			$ybucket = intval(floor($x[2] / $bucketsize));
	
			//print($xbucket.' '.$ybucket."\n");
	
			if (! array_key_exists($xbucket, $buckets))
			{
				$buckets[$xbucket] = array();
				$kdts[$xbucket] = array();
			}
	
			if (! array_key_exists($ybucket, $buckets[$xbucket]))
			{
				$buckets[$xbucket][$ybucket] = array();
				$kdts[$xbucket][$ybucket] = KdTree::make(2);
			}
	
			$cities[$x[0]] = array($x[1], $x[2], FALSE);
			$buckets[$xbucket][$ybucket][] = array($x[1], $x[2], $x[0], FALSE);
			$kdts[$xbucket][$ybucket]->add(array($x[1], $x[2]), $x[0]);
		}

		file_put_contents($fn, serialize(array($cities, $buckets, $kdts)));
	}

	function distxy($x1, $y1, $x2, $y2)
	{
		return sqrt((($x2 - $x1) * ($x2 - $x1)) + (($y2 - $y1) * ($y2 - $y1)));
	}

	function dist($id1, $id2)
	{
		global $cities;

		return distxy($cities[$id1][0], $cities[$id1][1], $cities[$id2][0], $cities[$id2][1]);
	}

	/*
	$total = 0;

	for ($i = 0; $i != $bnum; ++$i)
	{
		for ($j = 0; $j != $bnum; ++$j)
		{
			if (count($buckets[$i][$j]) < 8)
			{
				print('Bucket '.$i.' '.$j.' has too few nodes.'."\n");
			}

			$total += count($buckets[$i][$j]);
		}
	}

	print($total."\n");
	*/

	$city0 = NULL;

	$edges = array();
	$edgl = array(0 => array(), 1 => array());

	$paths = array();

	function fcp($i, $j)
	{
		global $city0, $cities, $buckets, $kdts, $edges, $edgl, $paths, $showPaths, $k0, $kc, $opt, $citycheck;

		if (! array_key_exists($i, $paths))
		{
			$paths[$i] = array();
		}

		if (! array_key_exists($j, $paths[$i]))
		{
			$paths[$i][$j] = array();
		}

		$pathd0 = 0;
		$path0 = NULL;
	
		//$edges = array();
	
		for ($jjj = 0; $jjj < 2; ++$jjj)
		{
			$path = array();
	
			$buedges = $edges;
			$buedgl = $edgl;

			$citynum = rand(0, count($buckets[$i][$j]) - 1);

			$city = $buckets[$i][$j][$citynum][2];

			$stpt = $city;
		
			$cities[$city][2] = TRUE;
		
			$path[] = $city;
		
			$dist = 0;
		
			//print('0 0 '.$city.' 0'."\n");
		
			for ($cnt = 1; $cnt < count($buckets[$i][$j]); ++$cnt)
			{
				$found = FALSE;
		
				$k = $k0 - $jjj;
		
				//print('Cur: '.$city."\n");

				$city0 = $city;
				$nns = $kdts[$i][$j]->findKnn(array($cities[$city][0], $cities[$city][1]), $k);
	
				$p = 0;

				foreach ($nns as $n => $nn)
				{
					//print($nn[0].' '.$nn[2]."\n");

					$p += (1 / ($nn[0] * $nn[0]));
					$nns[$n][3] = $p;
				}

				$r = (mt_rand() / mt_getrandmax()) * $p;

				foreach ($nns as $n => $nn)
				{
					if ($nn[3] < $r)
					{
						continue;
					}

					$cities[$nn[2]][2] = TRUE;
	
					$dist += $nn[0];
	
					$found = TRUE;
	
					$edges[$city.'-'.$nn[2]] = TRUE;
					$edges[$nn[2].'-'.$city] = TRUE;

					if (! array_key_exists($city, $edgl[$jjj]))
					{
						$edgl[$jjj][$city] = array();
					}

					if (! array_key_exists($nn[2], $edgl[$jjj]))
					{
						$edgl[$jjj][$nn[2]] = array();
					}

					$edgl[$jjj][$city][] = $nn[2];
					$edgl[$jjj][$nn[2]][] = $city;

					$city = $nn[2];
	
					$path[] = $city;
	
					if ($showPaths)
					{
						print($cnt.' '.$k.' '.$nn[2].' '.$dist."\n");
					}

					break;
				}

				if (! $found)
				{
					print('BLARGH!'."\n");

					break;

					die();
				}
	
				/**/
				//print('Optimizing... ['.$dist.'] ('.count($path).')'."\n");
	
				$iter = 1;
	
				//$optimize = TRUE;
				//$optimize = ($nn[0] >= 100) && (count($path) > 1024);
				$optimize = FALSE;

				//print(($optimize ? 'OKAY' : 'ABORT')."\n");
	
				$city0 = NULL;
	
				if ($optimize)
				{
					//$oldcities = $cities;

					$mindist = $dist;
					$minpath = $path;
					$minedges = $edges;
					$minedgl = $edgl;

					//for ($k = 0; $k < count($buckets[$i][$j]); ++$k)
					//{
					//	$cities[$buckets[$i][$j][$k][2]][2] = FALSE;
					//}
				}
	
				//$optimizationSteps = $opt;
				$optimizationSteps = 1;
	
				$permut = array();
	
				while ($optimize)
				{
					$nst = $path[count($path) - 2];
					$nen = $path[count($path) - 1];
					$ndst = dist($nst, $nen);
	
					print('O: '.$nst.' '.$nen.' '.$ndst."\n");
	
					$citycheck = FALSE;
					$nns = $kdts[$i][$j]->findKnn(array($cities[$nen][0], $cities[$nen][1]), $kc);
					$citycheck = TRUE;
	
					$bestgain = 0;
					$bestnn1 = NULL;
					$bestnn2 = NULL;
					//$besthash = NULL;
	
					//$pth = implode(',', $path);
	
					foreach ($nns as $nn)
					{
						$nn1 = $nn[2];
	
						if (! array_key_exists($nn1, $edgl[$jjj]))
						{
							continue;
						}

						foreach ($edgl[$jjj][$nn1] as $nn2)
						{
							//$hash = sha1($pth.':'.$nn1.':'.$nn2);
	
							//if (isset($edges[$nn1.'-'.$nen]) || isset($edges[$nn2.'-'.$nen]) || array_key_exists($hash, $permut))
							if (isset($edges[$nn1.'-'.$nen]) || isset($edges[$nn2.'-'.$nen]))
							{
								continue;
							}
	
							$gain = -(dist($nn1, $nen) + dist($nen, $nn2) - dist($nst, $nen) - dist($nn1, $nn2));
	
							//print('Gain: '.$gain.' ('.$nn1.', '.$nn2.')'."\n");
	
							if ($gain > $bestgain)
							{
								$bestgain = $gain;
								$bestnn1 = $nn1;
								$bestnn2 = $nn2;
								//$besthash = $hash;
							}
						}
					}
	
					print('Best gain: '.$bestgain.' ('.$bestnn1.', '.$bestnn2.')'."\n");
	
					if (is_null($bestnn1))
					{
						break;
					}
					else
					{
						$nn1ix = array_search($bestnn1, $path);
						$nn2ix = array_search($bestnn2, $path);
	
						$nn1aix = min($nn1ix, $nn2ix);
						$nn2aix = max($nn1ix, $nn2ix);
	
						if (($nn1aix + 1) !== ($nn2aix))
						{
							print('Invalid edge indices: '.$nn1aix.', '.$nn2aix.' ('.$bestnn1.', '.$bestnn2.')'."\n");
	
							die();
						}
	
						//$permut[$besthash] = TRUE;
	
						unset($edges[$nst.'-'.$nen]);
						unset($edges[$nen.'-'.$nst]);
						unset($edges[$bestnn1.'-'.$bestnn2]);
						unset($edges[$bestnn2.'-'.$bestnn1]);
						$edges[$bestnn1.'-'.$nen] = TRUE;
						$edges[$nen.'-'.$bestnn1] = TRUE;
						$edges[$bestnn2.'-'.$nen] = TRUE;
						$edges[$nen.'-'.$bestnn2] = TRUE;
	
						$edgl[$jjj][$nen] = array_filter($edgl[$jjj][$nen], function($x) use($nst) { return $x !== $nst; });
						$edgl[$jjj][$nst] = array_filter($edgl[$jjj][$nst], function($x) use($nen) { return $x !== $nen; });
						$edgl[$jjj][$bestnn1] = array_filter($edgl[$jjj][$bestnn1], function($x) use($bestnn2) { return $x !== $bestnn2; });
						$edgl[$jjj][$bestnn2] = array_filter($edgl[$jjj][$bestnn2], function($x) use($bestnn1) { return $x !== $bestnn1; });
						$edgl[$jjj][$nen][] = $bestnn1;
						$edgl[$jjj][$nen][] = $bestnn2;
						$edgl[$jjj][$bestnn1][] = $nen;
						$edgl[$jjj][$bestnn2][] = $nen;
	
						$path = array_merge(array_slice($path, 0, $nn1aix + 1), array($nen), array_slice($path, $nn2aix, count($path) - $nn1aix - 2));
	
						$dist -= $bestgain;
					}
	
					if ($dist < $mindist)
					{
						$mindist = $dist;
						$minpath = $path;
						$minedges = $edges;
						$minedgl = $edgl;
					}
	
					print($iter.': ['.$dist.'] ('.count($path).')'."\n");
	
					++$iter;
	
					//$optimize = FALSE;
	
					--$optimizationSteps;
	
					if ($optimizationSteps === 0)
					{
						break;
					}
				}
	
				if ($optimize)
				{
					//$cities = $oldcities;

					$dist = $mindist;
					$path = $minpath;
					$edges = $minedges;
					$edgl = $minedgl;
				}

				$city = $path[count($path) - 1];
				/**/
			}
		
			if (! $found)
			{
				$edges = $buedges;
				$edgl = $buedgl;

				for ($k = 0; $k < count($buckets[$i][$j]); ++$k)
				{
					$cities[$buckets[$i][$j][$k][2]][2] = FALSE;
				}

				--$jjj;

				continue;
			}

			/**/
			//print('Optimizing... ['.$dist.'] {'.implode(', ', $path).'}'."\n");
			print('Optimizing... ['.$dist.']'."\n");

			$iter = 1;

			$optimize = TRUE;

			$city0 = NULL;

			for ($k = 0; $k < count($buckets[$i][$j]); ++$k)
			{
				$cities[$buckets[$i][$j][$k][2]][2] = FALSE;
			}

			$mindist = $dist;
			$minpath = $path;
			$minedges = $edges;
			$minedgl = $edgl;

			$optimizationSteps = $opt;

			$permut = array();

			while ($optimize)
			{
				$nst = $path[count($path) - 2];
				$nen = $path[count($path) - 1];
				$ndst = dist($nst, $nen);

				//print($nst.' '.$nen.' '.$ndst."\n");

				$nns = $kdts[$i][$j]->findKnn(array($cities[$nen][0], $cities[$nen][1]), $kc);

				$bestgain = -100500;
				$bestnn1 = NULL;
				$bestnn2 = NULL;
				$besthash = NULL;

				$pth = implode(',', $path);

				foreach ($nns as $nn)
				{
					$nn1 = $nn[2];

					foreach ($edgl[$jjj][$nn1] as $nn2)
					{
						$hash = sha1($pth.':'.$nn1.':'.$nn2);

						if (isset($edges[$nn1.'-'.$nen]) || isset($edges[$nn2.'-'.$nen]) || array_key_exists($hash, $permut))
						{
							continue;
						}

						$gain = -(dist($nn1, $nen) + dist($nen, $nn2) - dist($nst, $nen) - dist($nn1, $nn2));

						//print('Gain: '.$gain.' ('.$nn1.', '.$nn2.')'."\n");

						if ($gain > $bestgain)
						{
							$bestgain = $gain;
							$bestnn1 = $nn1;
							$bestnn2 = $nn2;
							$besthash = $hash;
						}
					}
				}

				//print('Best gain: '.$bestgain.' ('.$bestnn1.', '.$bestnn2.')'."\n");

				if (is_null($bestnn1))
				{
					break;
				}
				else
				{
					$nn1ix = array_search($bestnn1, $path);
					$nn2ix = array_search($bestnn2, $path);

					$nn1aix = min($nn1ix, $nn2ix);
					$nn2aix = max($nn1ix, $nn2ix);

					if (($nn1aix + 1) !== ($nn2aix))
					{
						print('Invalid edge indices: '.$nn1aix.', '.$nn2aix.' ('.$bestnn1.', '.$bestnn2.')'."\n");

						die();
					}

					$permut[$besthash] = TRUE;

					unset($edges[$nst.'-'.$nen]);
					unset($edges[$nen.'-'.$nst]);
					unset($edges[$bestnn1.'-'.$bestnn2]);
					unset($edges[$bestnn2.'-'.$bestnn1]);
					$edges[$bestnn1.'-'.$nen] = TRUE;
					$edges[$nen.'-'.$bestnn1] = TRUE;
					$edges[$bestnn2.'-'.$nen] = TRUE;
					$edges[$nen.'-'.$bestnn2] = TRUE;

					$edgl[$jjj][$nen] = array_filter($edgl[$jjj][$nen], function($x) use($nst) { return $x !== $nst; });
					$edgl[$jjj][$nst] = array_filter($edgl[$jjj][$nst], function($x) use($nen) { return $x !== $nen; });
					$edgl[$jjj][$bestnn1] = array_filter($edgl[$jjj][$bestnn1], function($x) use($bestnn2) { return $x !== $bestnn2; });
					$edgl[$jjj][$bestnn2] = array_filter($edgl[$jjj][$bestnn2], function($x) use($bestnn1) { return $x !== $bestnn1; });
					$edgl[$jjj][$nen][] = $bestnn1;
					$edgl[$jjj][$nen][] = $bestnn2;
					$edgl[$jjj][$bestnn1][] = $nen;
					$edgl[$jjj][$bestnn2][] = $nen;

					$path = array_merge(array_slice($path, 0, $nn1aix + 1), array($nen), array_slice($path, $nn2aix, count($path) - $nn1aix - 2));

					$dist -= $bestgain;
				}

				if ($dist < $mindist)
				{
					$mindist = $dist;
					$minpath = $path;
					$minedges = $edges;
					$minedgl = $edgl;
				}

				//print($iter.': ['.$dist.'] ('.count($path).')'."\n");
				print($iter.': ['.$dist.'] Gain: ('.$bestgain.') ndst: ('.$ndst.')'."\n");

				++$iter;

				//$optimize = FALSE;

				--$optimizationSteps;

				if ($optimizationSteps === 0)
				{
					break;
				}
			}

			$dist = $mindist;
			$path = $minpath;
			$edges = $minedges;
			$edgl = $minedgl;
			/**/

			$stpt = $path[0];
			$endpt = $path[count($path) - 1];

			if (is_null($path0))
			{
				$path0 = $path;
				$pathd0 = $dist;
			}
	
			$paths[$i][$j][] = array($dist, $stpt, $endpt, $path);

			/*
			for ($k = 0; $k < count($buckets[$i][$j]); ++$k)
			{
				$cities[$buckets[$i][$j][$k][2]][2] = FALSE;
			}
			*/
		}
	
		print($pathd0.' '.$dist."\n");

		//die();
	}

	for ($i = 0; $i != $bnum; ++$i)
	{
		for ($j = 0; $j != $bnum; ++$j)
		{
			$cleanedges = $edges;

			fcp($i, $j);

			for ($k = 0; $k != ($runs - 1); ++$k)
			{
				$pp = $paths[$i][$j];
				$bestedges = $edges;

				$edges = $cleanedges;
				fcp($i, $j);

				if (($pp[0][0] + $pp[1][0]) < ($paths[$i][$j][0][0] + $paths[$i][$j][1][0]))
				{
					$paths[$i][$j] = $pp;
					$edges = $bestedges;
				}
			}

			//print($paths[$i][$j][0][0]."\n");
			//print($paths[$i][$j][1][0]."\n");
			//print_r($paths[$i][$j][0][3]);
		}
	}

	$path1 = array();
	$path2 = array();
	$dd1 = 0;
	$dd2 = 0;
	$lasti = 0;
	$lastj = 0;

	$needNotConnect = TRUE;

	for ($i = 0; $i != $bnum; ++$i)
	{
		for ($j = 0; $j != $bnum; ++$j)
		{
			$swap = FALSE;
			$p1 = $paths[$i][(($i % 2) == 0) ? $j : ($bnum - $j - 1)][($i + ((($i % 2) == 0) ? $j : ($bnum - $j - 1))) % 2];
			//$p1 = $paths[$i][(($i % 2) == 0) ? $j : ($bnum - $j - 1)][0];
			//print_r($p1);
			//if ((! $needNotConnect) && isset($edges[$cep1.'-'.$p1[1]]))
			if ((! $needNotConnect) && (isset($edges[$cep1.'-'.$p1[1]]) ||
					((! isset($edges[$cep1.'-'.$p1[2]])) && (dist($cep1, $p1[2]) < dist($cep1, $p1[1])))))
			{
				print($cep1.'-'.$p1[1].' already taken (1).'."\n");
				$swap = TRUE;
			}
			if ($swap)
			{
				$tmp = $p1[1];
				$p1[1] = $p1[2];
				$p1[2] = $tmp;
				$p1[3] = array_reverse($p1[3]);

				if (isset($edges[$cep1.'-'.$p1[1]]))
				{
					print('Alarm! Swapping didn\'t help. '.$cep1.'-'.$p1[1].' present as well. (1)'."\n");
					die();
				}
			}
			$path1 = array_merge($path1, $p1[3]);
			$dd1 += ($needNotConnect ? 0 : dist($cep1, $p1[1])) + $p1[0];

			if (isset($cep1))
			{
				$edges[$cep1.'-'.$p1[1]] = TRUE;
				$edges[$p1[1].'-'.$cep1] = TRUE;
			}

			$cep1 = $p1[2];

			$swap = FALSE;
			$p2 = $paths[(($i % 2) == 0) ? $j : ($bnum - $j - 1)][$i][($i + ((($i % 2) == 0) ? $j : ($bnum - $j - 1)) + 1) % 2];
			//$p2 = $paths[(($i % 2) == 0) ? $j : ($bnum - $j - 1)][$i][1];
			//print_r($p2);
			//if ((! $needNotConnect) && isset($edges[$cep2.'-'.$p2[1]]))
			if ((! $needNotConnect) && (isset($edges[$cep2.'-'.$p2[1]]) ||
					((! isset($edges[$cep2.'-'.$p2[2]])) && (dist($cep2, $p2[2]) < dist($cep2, $p2[1])))))
			{
				print($cep2.'-'.$p2[1].' already taken (2).'."\n");
				$swap = TRUE;
			}
			if ($swap)
			{
				$tmp = $p2[1];
				$p2[1] = $p2[2];
				$p2[2] = $tmp;
				$p2[3] = array_reverse($p2[3]);

				if (isset($edges[$cep2.'-'.$p2[1]]))
				{
					print('Alarm! Swapping didn\'t help. '.$cep1.'-'.$p1[1].' present as well. (2)'."\n");
					die();
				}
			}
			$path2 = array_merge($path2, $p2[3]);
			$dd2 += ($needNotConnect ? 0 : dist($cep2, $p2[1])) + $p2[0];

			if (isset($cep2))
			{
				$edges[$cep2.'-'.$p2[1]] = TRUE;
				$edges[$p2[1].'-'.$cep2] = TRUE;
			}

			$cep2 = $p2[2];

			print($i.' '.$j.' '.$cep1.' '.$cep2.' '.$dd1.' '.$dd2."\n");

			$needNotConnect = FALSE;

			$lasti = $i;
			$lastj = $j;
		}
	}

	print($dd1.' '.$dd2."\n");

	file_put_contents(__DIR__.'/sol.csv', 'path1,path2'."\n".implode("\n", array_map(function ($x, $y) { return $x.','.$y; }, $path1, $path2)));

?>
