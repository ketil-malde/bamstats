-- | BAM - extract information from BAM files

{-# Language DeriveDataTypeable #-}

module Main where
import Bio.SamTools.Classify
import Bio.SamTools.Bam
import System.Console.CmdArgs

data Options = Stats { numrd :: Maybe Int, inputs :: [FilePath] }
             | Hist { numrd :: Maybe Int, bins, maxdist :: Int, inputs :: [FilePath] }
             | Quants { numrd :: Maybe Int, mindist :: Int, delta :: Double, inputs :: [FilePath] }
             | Dump { numrd :: Maybe Int, inputs :: [FilePath] }
             deriving (Data,Typeable,Read,Show)

clfy, hst, quants, dump :: Options
clfy = Stats { numrd = Nothing &= help "max number of reads (default: all)"
                , inputs = []     &= args &= typ "BAM file(s)"
                } &= help "Calculate statistics on insert sizes"
hst  = Hist { numrd = Nothing &= help "max number of reads to include"
            , bins  = 20      &= help "number of bins" &= typ "INT"
            , maxdist = 1000  &= help "max insert size to count" &= typ "INT"
            , inputs = []     &= args &= typ "BAM file(s)"
            } &= help "Collect a histogram of insert sizes"
quants = Quants {
              numrd = Nothing &= help "max number of reads to include"
--            , bins  =  200    &= help "number of bins" &= typ "INT"
            , mindist = 100   &= help "minimum insert size to count" &= typ "INT"
            , delta   = 1.01  &= help "relative size of next bin" &= typ "FLOAT"                        
            , inputs = []     &= args &= typ "BAM file(s)"                
            } &= help "output approximate quantiles instead of the full histogram"
dump = Dump { numrd = Just 100 &= help "max number of reads (default: 100)" 
            , inputs = []     &= args &= typ "BAM file(s)"            
            } &= help "Dump alignments in the different classes."       

histgen :: Options -> Hist
histgen o = H 0 [(b,0) | b <- enumFromThenTo size (2*size) (maxdist o)]
  where size = (maxdist o `div` bins o)

quantgen :: Options -> Hist
quantgen o = H 0 [(round b,0) | b <- scanl1 (+) $ iterate (* delta o) (fromIntegral $ mindist o)]

main :: IO ()
main = do
  o <- cmdArgs $ modes [clfy,hst,quants,dump] 
       &= help "Extract information from BAM files" 
       &= program "bam" &= summary "bam v0.0, ©2012 Ketil Malde"
  let geninp f = (case numrd o of Just x -> take x; Nothing -> id) `fmap` readBams f
      genout = putStrLn . case o of Stats {} -> display . classify (CS 0 0 0 0 0)
                                    Hist {}  -> display . classify (histgen o)
                                    Quants {} -> showQuants . classify (quantgen o)
                                    Dump {}  -> display . classify (Bams [])
  mapM_ (\f -> genout =<< geninp f) $ inputs o