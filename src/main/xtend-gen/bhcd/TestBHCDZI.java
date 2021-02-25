package bhcd;

import bhcd.BHCD;
import bhcd.BasicDendrogram;
import bhcd.GraphAdjMatrix;
import blang.types.StaticUtils;
import blang.types.internals.RealScalar;
import blang.validation.ExactInvarianceTest;
import blang.validation.Instance;
import java.util.function.Function;
import org.eclipse.xtext.xbase.lib.ObjectExtensions;
import org.eclipse.xtext.xbase.lib.Procedures.Procedure1;
import org.junit.Test;

@SuppressWarnings("all")
public class TestBHCDZI {
  @Test
  public void invariance() {
    final GraphAdjMatrix.MutableGraphAdjMatrix graph = new GraphAdjMatrix.MutableGraphAdjMatrix(20);
    final RealScalar alpha = StaticUtils.latentReal();
    final RealScalar beta = StaticUtils.latentReal();
    final RealScalar theta = StaticUtils.latentReal();
    final BHCD model = new BHCD.Builder().setGraph(graph).setDend(BasicDendrogram.initDend(graph)).setAlpha(alpha).setBeta(beta).setTheta(theta).build();
    final Function<BHCD, Double> _function = (BHCD it) -> {
      double _xifexpression = (double) 0;
      int _numEdges = graph.getNumEdges();
      int _numVertices = graph.getNumVertices();
      boolean _greaterThan = (_numEdges > _numVertices);
      if (_greaterThan) {
        _xifexpression = 1.0;
      } else {
        _xifexpression = 0.0;
      }
      return Double.valueOf(_xifexpression);
    };
    final Function<BHCD, Double> _function_1 = (BHCD it) -> {
      return Double.valueOf(it.getAlpha().doubleValue());
    };
    final Function<BHCD, Double> _function_2 = (BHCD it) -> {
      return Double.valueOf(it.getBeta().doubleValue());
    };
    final Function<BHCD, Double> _function_3 = (BHCD it) -> {
      return Double.valueOf(it.getTheta().doubleValue());
    };
    final Function<BHCD, Double> _function_4 = (BHCD it) -> {
      int _numEdges = graph.getNumEdges();
      int _numVertices = graph.getNumVertices();
      int _numVertices_1 = graph.getNumVertices();
      int _minus = (_numVertices_1 - 1);
      int _multiply = (_numVertices * _minus);
      int _divide = (_multiply / 2);
      return Double.valueOf(Integer.valueOf((_numEdges / _divide)).doubleValue());
    };
    final Instance<BHCD> instance = new Instance<BHCD>(model, _function, _function_1, _function_2, _function_3, _function_4);
    ExactInvarianceTest _exactInvarianceTest = new ExactInvarianceTest();
    final Procedure1<ExactInvarianceTest> _function_5 = (ExactInvarianceTest it) -> {
      it.add(instance);
    };
    final ExactInvarianceTest test = ObjectExtensions.<ExactInvarianceTest>operator_doubleArrow(_exactInvarianceTest, _function_5);
    test.nPosteriorSamplesPerIndep = 20;
    test.nIndependentSamples = 20_000;
    test.check();
  }
}
