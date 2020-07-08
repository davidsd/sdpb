#include <El.hpp>

using namespace std::literals;

std::vector<El::BigFloat>
get_new_points(const std::vector<El::BigFloat> &points)
{
  std::vector<std::string> p(
    {"0.0"s, "100.0"s,
     "70.71097310187551316924294967198658288202601723472274987600171978175689041069736090738990862176715189781008228550345393015742895478251094421336161102868679079043214338263992547395032484562279010968063798428727009843"s,
     "50.00062625468106944019216182643156474379677118076837100467537272049853550244517869632018202076432364570277375378184317259970870495580695072169212498074832410294199146110624256140190310179358112574339162711420578657"s,
     "35.35637397962831347683863823502104974040889756157254545482681650663951276898283472697374576567768192516309111325368095659517475221328154328379409112992216968073938576240924428432260294664860391313909042244905654799"s,
     "25.00157314902451015465465402635618190797057537266955746650724265360108030824253482783000938914219283163388934810954078777881444673809186370835417669438108874736147051659167910815046145405458225799344593551677712334"s,
     "17.67998309320581741800427653299058748336717297759677327965699280422733563295723146765404049145597552131172173975293390990260802256438167680318260613219067359568857670125749362496290002305028338351107736873445511594"s,
     "12.50336670187752790870980618019030772510207089325191592290938910695487695214286329240325580643832288362256162437532732518884403505092153429414786923848352737604954521476423251438585769002046109656758198038610234114"s,
     "8.84375363277898848681983744549812298616518854387200491066307577895636542142687717372969887421684891050779517691560391155510132658758345295916048258853481735284930389850459504969098552503587172451073321426601498035"s});
  std::vector<El::BigFloat> result;
  if(points.size() < p.size())
    {
      result.emplace_back(p[points.size()]);
    }
  return result;
}