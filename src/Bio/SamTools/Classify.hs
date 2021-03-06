
-- | Stats reads in a BAM-file as innie outie rightie leftie,
--   with stats for distances

module Bio.SamTools.Classify where

import Prelude hiding (sum)
import Bio.SamTools.Bam
import Data.List (foldl', intercalate)
import Data.Maybe (isNothing, fromJust)
import Text.Printf (printf)

data Stats a = Class { total, unmapped, orphans, splits, secondary :: !Int, innies, outies, lefties, righties :: !a } deriving Show

-- --------------------------------------------------
-- The generic collection framework
-- --------------------------------------------------

-- | Pretty-print a 'Stats' value.
display :: Insertable x => Stats x -> String
display c = unlines $ map (intercalate "\t")
  [  "#Alignment" : dispheader (innies c)
  ,  "innies   " : disp1 t (innies c)
  ,  "outies   " : disp1 t (outies c)
  ,  "lefties  " : disp1 t (lefties c)
  ,  "righties " : disp1 t (righties c)
  , summarize c t]
  where t = total c - secondary c

summarize :: Stats a -> Int -> [String]
summarize c t = [printf "\nTotal reads:  %7d\n" t
    ++" unmapped:    "++show1 unmapped++"\n"
    ++" orphans:     "++show1 orphans++"\n"
    ++" split pairs: "++show2 splits++"\n"
    ++" secondary:   "++show1 secondary]
   where show1 f = printf "%7d" (f c)++percent (f c)
         show2 f = printf "%7d" (f c`div`2)++percent (f c)
         percent x = printf " (%.1f%%)" (100.0*fromIntegral x/fromIntegral t::Double)

showQuants :: Stats Hist -> String
showQuants cs = unlines $ map (intercalate "\t")
  [  "#Alignment" : hdr (innies cs)
  ,  "innies   " : quants t (innies cs)
  ,  "outies   " : quants t (outies cs)
  ,  "lefties  " : quants t (lefties cs)
  ,  "righties " : quants t (righties cs)
  , summarize cs t]
  where t = total cs - secondary cs
        percentiles = [0.05,0.25,0.5,0.75,0.95] :: [Double]
        quants tot h = printf "%14d" (hcount h)
              : printf "%5.2f%%" (200*fromIntegral (hcount h)/fromIntegral tot::Double)
              : if hcount h == 0 then map (const "    N/A") percentiles
                else go (map (round . (*fromIntegral (hcount h))) percentiles) 0 (0,0) (buckets h)
        hdr _h = "         count":"  prop":map (printf "%6.0f%%" . (*100)) percentiles
        go :: [Int] -> Int -> (Int,Int) -> [(Int,Int)] -> [String]
        go (0:fs) sum prev inp = "      0" : go fs sum prev inp
        go (f:fs) sum prev ((b,c):rest) = 
          let next = sum+c in 
          if next >= f
          then let b' = b - round (fromIntegral (next-f)/fromIntegral c*fromIntegral (b-fst prev)::Double)
               in printf "%7d" b' : (go fs sum prev ((b,c):rest))
          else     go (f:fs) next (b,c) rest
        go [] _ _ _ = []
        go _ _ _ [] = error "ran out of reads!?"

genplot :: String -> Stats Hist -> String
genplot cmds h = preamble ++ cmds ++ "\n" ++ plot
  where preamble = unlines ["set style data boxes"
                           ,"set style fill solid border 0"
                           ,"set xlabel 'insert length'"
                           ,"set ylabel 'read count'"]
        plot = "plot '-' ti 'innies', '-' ti 'outies', '-' ti 'lefties', '-' ti 'righties'\n"
                  ++concat [ plot1 $ buckets $ innies h
                           , plot1 $ map (\(x,y)->(-x,-y)) $ buckets $ outies h
                           , plot1 $ map (\(x,y)->(-x, y)) $ buckets $ lefties h
                           , plot1 $ map (\(x,y)->( x,-y)) $ buckets $ righties h]
        plot1 = unlines . go 0
        go prev ((b,c):rest@(_:_)) = (show ((prev+b) `div` 2) ++ " " ++ show c) : go (b+1) rest
        go _ [_] = ["e"] -- ignore last bucket, which is a catch-all
        go _ [] = error "no data"

-- | Something that collects information from Bam1s, similar
--   to a fold, and can display it as lines.
class Insertable x where
  insert :: Bam1 -> x -> x         -- collect from a Bam1
  disp1  :: Int  -> x -> [String]  -- the int is the total, for statistics
  dispheader :: x -> [String]      -- display the appropriate header

-- | Extract info from alignments
classify :: Insertable x => x -> [Bam1] -> Stats x
classify def = foldl' (class1 . bump) (Class 0 0 0 0 0 def def def def)

-- | Update count
bump :: Stats a -> Stats a
bump s = s { total = total s + 1 }

-- | Update data structure with a single alignment
class1 :: Insertable x => Stats x -> Bam1 -> Stats x
class1 c0 b
  | isSecondary b = add_secondary b c0
  | isNothing (insertSize b) = add_unmapped b c0
  | isOpposite b =
      (if firstUpstream b then add_innie else add_outie) b c0
  | otherwise =
      (if firstUpstream b then add_rightie else add_leftie) b c0

add_innie, add_outie, add_rightie, add_leftie, add_unmapped, add_secondary
  :: Insertable x => Bam1 -> Stats x -> Stats x
add_innie b c = c { innies = insert b (innies c) }
add_outie b c = c { outies = insert b (outies c) }
add_rightie b c = c { righties = insert b (righties c) }
add_leftie b c = c { lefties = insert b (lefties c) }

add_secondary _ c = c { secondary = secondary c + 1 }
add_unmapped b c  
  | isUnmap b                    = c { unmapped = unmapped c + 1 }
  | isMateUnmap b                = c { orphans = orphans c + 1 }
  -- both mateTgtID and tgtId should be Just something  
  | mateTargetID b /= targetID b = c { splits = splits c + 1 }
  | otherwise = c

isOpposite, firstUpstream, isBefore :: Bam1 -> Bool
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
statistics (Class t ump orp sp sec i o l r) = (Class t ump orp sp sec (mkstats i) (mkstats o) (mkstats l) (mkstats r)) -- todo: Functor instance?

instance Insertable ClassStats where
  insert b (CS c s s2 s3 s4) = CS (c+1) (s+d) (s2+d^(2::Int)) (s3+d^(3::Int)) (s4+d^(4::Int)) 
    where d = maybe (error ("no insert size?\n"++show b)) fromIntegral $ insertSize b
  dispheader _ = ["         count","   prop","   mean","  stdev","   skew","   kurt"]
  disp1 tot cs = printf "%14d" (ccount cs) 
               : printf "%5.2f%%" (200*fromIntegral (ccount cs)/fromIntegral tot::Double)
               : map (printf "%7.1f") [mean s, stdev s, skew s, kurt s]
    where 
      s = mkstats cs  

-- --------------------------------------------------
-- Hist for accumulating histograms of sizes
-- --------------------------------------------------

data Hist = H { hcount :: !Int, buckets :: ![(Int,Int)] } deriving Show

instance Insertable Hist where
  insert b h = H (hcount h + 1) (go (fromIntegral $ fromJust $ insertSize b) (buckets h))
    where go x ((b1,v):bs@(_:_)) =
            if x > b1 then (b1,v):go x bs else inc (b1,v):bs
          go _ [(b1,v)] = [inc (b1,v)]
          go _ [] = error "this never happens"
          inc (x,y) = let y' = y+1 in y' `seq` (x,y')
  disp1 tot h = printf "%14d" (hcount h)
              : printf "%5.2f%%" (200*fromIntegral (hcount h)/fromIntegral tot::Double)
              : map (printf "%7d" . snd) (buckets h)
  dispheader h = "         count":"  prop":(map (printf "%7d" . fst) $ init $ buckets h)++["    >"]

-- --------------------------------------------------
-- Just "Collect" all the Bams in different classes
-- --------------------------------------------------

data Collect = Bams { bams :: [Bam1] }

instance Insertable Collect where
  insert b (Bams bs) = Bams (b:bs)
  disp1 _ (Bams bs) = [unlines $ ("":map show bs)]
  dispheader = const ["foo"]
