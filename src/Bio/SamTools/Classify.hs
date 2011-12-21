
-- | Stats reads in a BAM-file as innie outie rightie leftie,
--   with stats for distances

module Bio.SamTools.Classify where

import Bio.SamTools.Bam
import Data.List (foldl', intercalate)
import Data.Maybe (isNothing, fromJust)
import Text.Printf (printf)

data Stats a = Class { total :: !Int, innies, outies, lefties, righties :: !a } deriving Show

-- --------------------------------------------------
-- The generic collection framework
-- --------------------------------------------------

-- | Pretty-print a 'Stats' value.
display :: Insertable x => Stats x -> String
display c = unlines $ map (intercalate "\t")
  [  "Alignment" : dispheader (innies c)
  ,  "innies   " : disp1 t (innies c)
  ,  "outies   " : disp1 t (outies c)
  ,  "lefties  " : disp1 t (lefties c)
  ,  "righties " : disp1 t (righties c)
  , ["Total reads: ", show t]]
  where t = total c

class Insertable x where
  insert :: Bam1 -> x -> x
  disp1  :: Int  -> x -> [String]
  dispheader :: x -> [String]
  cdef :: x  

-- | Extract info from alignments
classify :: Insertable x => x -> [Bam1] -> Stats x
classify def = foldl' (class1 . bump) (Class 0 def def def def)

-- | Update count
bump :: Stats a -> Stats a
bump s = s { total = total s + 1 }

-- | Update data structure with a single alignment
class1 :: Insertable x => Stats x -> Bam1 -> Stats x
class1 c0 b
  | isUnmapped b = c0
  | isOpposite b =
      (if firstUpstream b then add_innie else add_outie) b c0
  | otherwise =
      (if firstUpstream b then add_rightie else add_leftie) b c0

add_innie, add_outie, add_rightie, add_leftie :: Insertable x => Bam1 -> Stats x -> Stats x
add_innie b c = c { innies = insert b (innies c) }
add_outie b c = c { outies = insert b (outies c) }
add_rightie b c = c { righties = insert b (righties c) }
add_leftie b c = c { lefties = insert b (lefties c) }

isUnmapped, isOpposite, firstUpstream, isBefore :: Bam1 -> Bool
isUnmapped = isNothing . insertSize
-- isUnmap b || isMateUnmap b || mateTargetID b /= targetID b -- not the same
-- insertSize is Nothing if negative (downstream read)
isOpposite b = isReverse b /= isMateReverse b 
firstUpstream b = isReverse b /= isBefore b
isBefore b = position b < matePosition b

-- --------------------------------------------------
-- ClassStats for calculating mean, stdev, etc etc
-- --------------------------------------------------

data ClassStats = CS { ccount :: !Int, xsum, x2sum, x3sum, x4sum :: !Double } 
                deriving Show

data Statistics = S { ncount, mean, stdev, skew, kurt :: Double }

mkstats :: ClassStats -> Statistics
mkstats cs = S n m s sk kt
  where
    n = fromIntegral (ccount cs)
    x = xsum cs
    x2 = x2sum cs
    x3 = x3sum cs
    x4 = x4sum cs
    m = x/n
    m2 = m*m
    m3 = m*m2
    s = sqrt ((x2-m2*n)/(n-1))
    sk = (x3 - 3*m*x2 + 2*m3*n)/(s*s*s*n)
    kt = (x4 - 4*m*x3 + 6*m2*x2 - 4*m3*x + n*m*m3)/(s*s*s*s*n) - 3

statistics :: Stats ClassStats -> Stats Statistics
statistics (Class t i o l r) = (Class t (mkstats i) (mkstats o) (mkstats l) (mkstats r)) -- todo: Functor instance?

instance Insertable ClassStats where
  insert b (CS c s s2 s3 s4) = CS (c+1) (s+d) (s2+d^(2::Int)) (s3+d^(3::Int)) (s4+d^(4::Int)) 
    where d = maybe (error ("no insert size?\n"++show b)) fromIntegral $ insertSize b
  dispheader _ = ["         count","   prop","   mean","  stdev","   skew","   kurt"]
  cdef = CS 0 0 0 0 0
  disp1 tot cs = printf "%14d" (ccount cs) 
               : printf "%5.1f%%" (100*fromIntegral (ccount cs)/fromIntegral tot::Double)
               : map (printf "%7.1f") [mean s, stdev s, skew s, kurt s]
    where 
      s = mkstats cs  

-- --------------------------------------------------
-- Hist for accumulating histograms of sizes
-- --------------------------------------------------

data Hist = H { hcount :: !Int, buckets :: ![(Int,Int)] } deriving Show

instance Insertable Hist where
  insert b h = H (hcount h + 1) (go (fromIntegral $ fromJust $ insertSize b) (buckets h))
    where go x ((b1,v):bs@(_:_)) =   -- fixme: strictness!
            if x > b1 then (b1,v):go x bs else (b1,v+1):bs
          go x [(b1,v)] = [(b1,v+1)]
  disp1 tot h = printf "%14d" (hcount h)
              : printf "%5.1f%%" (100*fromIntegral (hcount h)/fromIntegral tot::Double)
              : map (printf "%7.1f" . (fromIntegral :: (Int->Double)) . snd) (buckets h)
  dispheader h = "         count":"  prop":(map (printf "%7d" . fst) $ init $ buckets h)++["    >"]
  cdef = undefined

-- --------------------------------------------------
-- Just "Collect" all the Bams in different classes
-- --------------------------------------------------

data Collect = Bams { bams :: [Bam1] }

instance Insertable Collect where
  insert b (Bams bs) = Bams (b:bs)
  disp1 _ (Bams bs) = [unlines $ ("":map show bs)]
  dispheader = const ["foo"]
  cdef = Bams []
